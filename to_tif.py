'''
Utility script to gdal_translate jpgs to tifs
'''

import re
import sys

from pathlib import Path

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


def translate():
    '''
    Translate all the jpgs in jpeg_root to jpeg-compressed tiffs in tif_root,
    maintaining folder structure.
    '''

    jpeg_root = r'c:\gis\projects\sanborn\marriott_source'
    tif_root = r'c:\gis\projects\sanborn\marriott_tif'

    year_regex = '[0-9]{4}'

    #: Set this option to prevent jpeg size error
    gdal.SetConfigOption('GDAL_ALLOW_LARGE_LIBJPEG_MEM_ALLOC', 'YES')

    # tif_path = Path(tif_root)
    # if not tif_path.exists():
    #     tif_path.mkdir()

    jpgs = sorted(Path(jpeg_root).glob('**/*.jpg'))

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

            print(destination_file)
            #: Translate .jpg to .tif with gdal.Translate
            creation_opts = ['tiled=yes', 'compress=jpeg', 'photometric=ycbcr',
                             'jpeg_quality=100']
            trans_opts = gdal.TranslateOptions(format='GTiff',
                                               creationOptions=creation_opts,
                                               callback=gdal_progress_callback)
            dataset = gdal.Translate(str(destination_file), 
                                     str(jpg), options=trans_opts)
            dataset = None

if __name__ == '__main__':
    translate()
