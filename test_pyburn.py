import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import unittest
import numpy as np
import os
from osgeo import gdal, osr, ogr  # type: ignore
from shapely.geometry import Polygon, box
import geopandas as gpd
import rasterio
from shapely.geometry import mapping
from rasterio.transform import from_origin
import fiona # type: ignore
import matplotlib.pyplot as plt
from pyburn import read_band_image, nbr, dnbr, array2raster, clip_raster, reclassify, reproject_shp_gdal

class TestPyburn(unittest.TestCase):

    def setUp(self):
        """
        Set up test environment by creating dummy raster and shapefile data.
        """
        # Create dummy raster
        self.raster_file = 'test_raster.tif'
        self.raster_data = np.array([[1, 2], [3, 4]], dtype=np.float32)
        self.transform = from_origin(0, 2, 1, 1)  # Top-left corner at (0, 2), pixel size = 1x1
        self.crs = 'EPSG:4326'

        with rasterio.open(
            self.raster_file,
            'w',
            driver='GTiff',
            height=self.raster_data.shape[0],
            width=self.raster_data.shape[1],
            count=1,
            dtype=self.raster_data.dtype,
            crs=self.crs,
            transform=self.transform,
        ) as dst:
            dst.write(self.raster_data, 1)

        # Create dummy shapefile
        self.shapefile = 'test_shapefile.shp'
        polygon = Polygon([(0.5, 1.5), (1.5, 1.5), (1.5, 0.5), (0.5, 0.5)])  # Small square
        gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs=self.crs)
        gdf.to_file(self.shapefile)

    def tearDown(self):
        """
        Clean up test environment by removing dummy files.
        """
        if hasattr(self, 'raster_dataset') and self.raster_dataset:
            self.raster_dataset.close()  # Ensure the raster dataset is closed
        if os.path.exists(self.raster_file):
            os.remove(self.raster_file)
        if os.path.exists(self.shapefile):
            for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg']:
                try:
                    os.remove(self.shapefile.replace('.shp', ext))
                except FileNotFoundError:
                    pass

    def test_read_band_image(self):
        band_path = self.raster_file  # Use the raster file created in setUp
        band_data, spatial_ref, geo_transform, target_prj = read_band_image(band_path)

        # Check if the band data is a numpy array
        self.assertIsInstance(band_data, np.ndarray)

        # Check if the metadata are returned as expected
        self.assertIsInstance(spatial_ref, str)  # Spatial reference should be a string
        self.assertIsInstance(geo_transform, tuple)  # GeoTransform should be a tuple
        self.assertIsInstance(target_prj, osr.SpatialReference)  # Target projection should be an osr.SpatialReference object

    def test_clip_raster(self):
        raster = np.array([[1, 2], [3, 4]], dtype=float)
        geoTransform = (0, 1, 0, 0, 0, -1)
        projection = "EPSG:4326"
        raster_path = 'test_raster.tif'
        shapefile_path = 'test_shapefile.shp'

        dataset, _ = array2raster(raster, geoTransform, projection, raster_path)
    
        # Create a small square polygon in shapefile matching one raster pixel
        schema = {'geometry': 'Polygon', 'properties': {}}
        with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', schema=schema, crs=projection) as shp:
            shp.write({'geometry': mapping(box(0, -1, 1, 0)), 'properties': {}})

        fire_boundary = gpd.read_file(shapefile_path)
        clipped, _, _, _ = clip_raster(raster_path, fire_boundary)
    
        # Ensure the clipped array dimensions match
        self.assertEqual(clipped.shape[1:], (1, 1))  # Ignore the band dimension if present
    
        dataset = None  # Ensure cleanup
        os.remove(raster_path)
        os.remove(shapefile_path)

    def test_array2raster(self):
        output_file = 'test_output.tif'
        array = np.array([[1, 2], [3, 4]])
        geoTransform = (0, 1, 0, 0, 0, -1)
        projection = "EPSG:4326"
    
        dataset, _ = array2raster(array, geoTransform, projection, output_file)
    
        self.assertTrue(os.path.exists(output_file))
    
        dataset = None  # Properly close the dataset
        os.remove(output_file)

    def test_reclassify(self):
        array = np.array([[0.05, 0.15], [0.3, 0.5]])  # Input array
        # Updated expected array to match function logic
        expected = np.array([[1, 2], [3, 4]])

        reclassified = reclassify(array).astype(int)

        np.testing.assert_array_equal(reclassified, expected)

    def test_nbr(self):
        """
        Test the NBR calculation.
        """
        band1 = np.array([[3, 6], [9, 12]])
        band2 = np.array([[1, 2], [3, 4]])
        result = nbr(band1, band2)

        # Expected NBR values
        expected = (band1 - band2) / (band1 + band2)
        np.testing.assert_almost_equal(result, expected)

    def test_dnbr(self):
        """
        Test the dNBR calculation.
        """
        nbr1 = np.array([[0.6, 0.8], [0.9, 1.0]])
        nbr2 = np.array([[0.5, 0.7], [0.8, 0.9]])
        result = dnbr(nbr1, nbr2)

        # Expected dNBR values
        expected = nbr1 - nbr2
        np.testing.assert_almost_equal(result, expected)

    def test_reproject_shp_gdal(self):
        """
        Test the reproject_shp_gdal function.
        """
        infile_shp = 'test_shapefile.shp'
        outfile_shp = 'test_reprojected_shapefile.shp'

        # Create a dummy shapefile
        polygon = Polygon([(0.5, 1.5), (1.5, 1.5), (1.5, 0.5), (0.5, 0.5)])  # Small square
        gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs='EPSG:4326')
        gdf.to_file(infile_shp)

        # Reproject the shapefile to a new CRS (e.g., EPSG:3857)
        targetprj = osr.SpatialReference()
        targetprj.ImportFromEPSG(3857)  # Reproject to EPSG:3857
        reproject_shp_gdal(infile_shp, outfile_shp, targetprj)

        # Check if the reprojected shapefile exists
        self.assertTrue(os.path.exists(outfile_shp))

        # Cleanup
        os.remove(infile_shp)
        for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg']:
            try:
                os.remove(outfile_shp.replace('.shp', ext))
            except FileNotFoundError:
                pass

if __name__ == '__main__':
    unittest.main()
