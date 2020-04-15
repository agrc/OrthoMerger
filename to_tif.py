'''
Utility script to gdal_translate jpgs to tifs
'''

import re
import sys

from pathlib import Path

import pandas as pd

from osgeo import gdal


#: GDAL callback method that seems to work, as per
#: https://gis.stackexchange.com/questions/237479/using-callback-with-python-gdal-rasterizelayer
def gdal_progress_callback(complete, message, unknown):
    '''
    Progress bar styled after the default GDAL progress bars. Uses specific
    signature to conform with GDAL core.
    '''
    #: 40 stops on our progress bar, so scale to 40
    done = int(40 * complete / 1)

    #: Build string: 0...10...20... - done.
    status = ''
    for i in range(0, done):
        if i % 4 == 0:
            status += str(int(i / 4 * 10))
        else:
            status += '.'
    if done == 40:
        status += '100 - done.\n'

    sys.stdout.write('\r{}'.format(status))
    sys.stdout.flush()
    return 1


def get_projection(file_path):

    file_handle = gdal.Open(file_path, gdal.GA_ReadOnly)
    projection = file_handle.GetProjection()
    file_handle = None

    return projection

def get_resolution(file_path):

    file_handle = gdal.Open(file_path, gdal.GA_ReadOnly)
    transform = file_handle.GetGeoTransform()
    raster_xwidth = transform[1]
    raster_yheight = transform[5]

    return (raster_xwidth, raster_yheight)


def translate():
    '''
    Translate all the jpgs in jpeg_root to jpeg-compressed tiffs in tif_root,
    maintaining folder structure.
    '''

    jpeg_root = r'c:\gis\projects\sanborn\marriott_source_test'
    tif_root = r'c:\gis\projects\sanborn\marriott_tif_test'

    year_regex = '[0-9]{4}'

    #: Set this option to prevent jpeg size error
    gdal.SetConfigOption('GDAL_ALLOW_LARGE_LIBJPEG_MEM_ALLOC', 'YES')
    gdal.UseExceptions()

    jpgs = sorted(Path(jpeg_root).glob('**/*.jpg'))

    counter = 1
    error_files = []
    projection_dict = {}
    projection_list = []
    resolution_dict = {}

    for jpg in jpgs:

        #: Make sure it exists in a year subdiretory
        if re.search(year_regex, jpg.parent.name):

            #: Make sure we have destination directory
            year = jpg.parent.name
            city = jpg.parents[1].name
            destination_dir = Path(tif_root, city, year)
            if not destination_dir.exists():
                destination_dir.mkdir(parents=True)

            #: Create destination file name
            tif_name = jpg.stem + '.tif'
            destination_file = destination_dir / tif_name

            print(f'{counter} of {len(jpgs)}: {destination_file}')

            #: Translate .jpg to .tif with gdal.Translate
            # creation_opts = ['tiled=yes', 'compress=jpeg', 'photometric=ycbcr',
            #                  'jpeg_quality=100']
            creation_opts = ['tiled=yes', 'compress=lzw']
            try:
                # trans_opts = gdal.TranslateOptions(format='GTiff',
                #                                 creationOptions=creation_opts,
                #                                 callback=gdal_progress_callback)
                # dataset = gdal.Translate(str(destination_file),
                #                         str(jpg), options=trans_opts)
                # dataset = None


                #: Manually set source projection based on text of source file's
                #: projection property.
                projection = get_projection(str(jpg))
                #: Default to utm 12n unless it's blank or utah central (all
                #: blank projections are currently utah central)
                source_crs = 'epsg:26912'
                if not projection or 'Utah Central (ftUS)' in projection:
                    source_crs = 'epsg:3566'

                destination_crs = 'epsg:26912'

                warp_opts = gdal.WarpOptions(srcSRS=source_crs,
                                             dstSRS=destination_crs,
                                             resampleAlg='cubic',
                                             srcNodata='255 255 255',
                                             dstNodata='256 256 256',
                                             format='GTiff',
                                             outputType=gdal.GDT_Int16,
                                             multithread=True,
                                             creationOptions=creation_opts,
                                             callback=gdal_progress_callback)
                dataset = gdal.Warp(str(destination_file), str(jpg), options=warp_opts)
                dataset = None  #: Releases file handle

                # projection_text = get_projection(str(jpg))
                # projection_dict[str(jpg)] = projection_text
                # if projection_text not in projection_list:
                #     projection_list.append(projection_text)

            except RuntimeError as e:
                error_files.append(e)
                print(e)

            counter += 1

    if error_files:
        print('Failed files:')
        for e in error_files:
            print(e)

    # projection_frame = pd.DataFrame.from_dict(projection_dict, orient='index')
    # csv_path = Path(tif_root, 'projections.csv')
    # projection_frame.to_csv(csv_path)
    # for p in projection_list:
    #     print(p)

if __name__ == '__main__':
    translate()
