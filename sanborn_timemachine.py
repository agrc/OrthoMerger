'''
Explore Sanborn data delivered in a directory structure of root/city/year
(the rasters should have been reprojected and translated to tif previously).
'''

import datetime

from pathlib import Path

import merge
import re

def time_travel(main_dir_path, output_path):
    '''
    Walk through main_dir_path and build a list of vintages (city+year). Feed
    vintages one at a time to merge.run() to build a single mosaiced tif of
    that vintage.
    '''

    start = datetime.datetime.now()
    problems = []
    problem_years = []

    year_regex = '\d{4}$'

    #: Build list of paths to each vintage in main_dir path
    #: (.../main_dir_path/city/year)
    vintages = []
    for city_dir in Path(main_dir_path).iterdir():
        for year_dir in city_dir.iterdir():
            #: Only add directories with four-number year name (ie, '1911')
            if re.match(year_regex, year_dir.name):
                vintages.append(year_dir)

    #: Process vintages
    for year_dir in vintages:
        year = year_dir.name
        city = year_dir.parent.name
        filename_root = f'{city}{year}'
        output_dir = Path(output_path, city)

        try:
            merge.run(year_dir, output_dir, filename_root, fishnet_size=10, cleanup=False, tile=True)
        except Exception as e:
            problems.append(e)
            problem_years.append(f'{city}{year}')

    if problem_years:
        print(problems)
        print(problem_years)

    end = datetime.datetime.now()
    print(f'\nTotal run time: {end-start}')


if __name__ == '__main__':
    source = r'C:\gis\Projects\Sanborn\marriott_tif'
    destination = r'F:\WasatchCo\sanborn2'
    time_travel(source, destination)
