# orthomerger

Merges overlapping orthorectified rasters (imagery, scanned maps, etc) into a single output file.

Currently, `orthomerger` is a workable, if clunky, mosaicing utility for rasters. `merge.py`'s `run()` is the main entry into the program. `sanborn_timemachine.py` is a project-specific example of how it can be run while maintaining a specific directory structure.

In brief, it slices all the source rasters into individual tiles conforming to a single fishnet grid and sorts them based on distance to the center of the source raster and number of nodata cells. The sorted tiles are combined into a single VRT, which is then translated into a single 8-bit, jpeg-compressed GeoTiff (with overviews).

`orthomerger` works best with rasters of similar size that have overlapping data areas, such as an aerial imagery flight line. It does not work well for source rasters with irregular collars, like USGS topos with a larger collar area on bottom than top.

`orthomerger` thinks in terms of "source rasters", "chunks", "tiles", and "mosaics"

**source raster:**
    The initial rasters that will be mosaiced together. They should all live in the same directory and have a proper nodata value set. If working with 8-bit aerial imagery, it may be necessary to convert the source rasters into a 16-bit file and assign a nodata value prior to any warping/reprojecting.

**chunk:**
    A single, specific area of the fishnet grid used to slice all of the source rasters into individual tiles. Each chunk has a row/col index. The output `_mosaic.shp` shapefile contains this fishnet grid, with coincident polygons for every source raster that covers a specific chunk.

**tile:**
    A slice of a specific source raster for a specific chunk. Each tile has information used to help sort the final collection of tiles for each chunk to determine which should be on top.

**mosaic:**
    The final output. Hopefully a gorgeously tiled version of all the source rasters.

After performing an intial mosaic, the user can edit the resulting `_mosaic.shp` shapefile to overwrite the tiling algorithm and force a specific tile to the top. To do so, add a `Y` to the `overwrite` field for the polygon representing the desired tile.
