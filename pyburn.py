import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rasterio
import glob
from rasterio.plot import show
from rasterio.mask import mask
from scipy import datasets
from shapely.geometry import mapping
from osgeo import gdal, osr, ogr #type: ignore
gdal.UseExceptions()
import geopandas as gpd
import math
from matplotlib.colors import Normalize  # Import Normalize

# Function definitions

def read_band_image(band_path):
    """
    Reads the band image and returns the data as an array and associated metadata.
    """
    try:
        img = gdal.Open(band_path)
        if img is None:
            raise FileNotFoundError(f"Failed to open {band_path}")
        data = np.array(img.GetRasterBand(1).ReadAsArray())
        spatialRef = img.GetProjection()
        geoTransform = img.GetGeoTransform()
        targetprj = osr.SpatialReference(wkt=img.GetProjection())
        return data, spatialRef, geoTransform, targetprj
    except Exception as e:
        raise RuntimeError(f"Error reading band image {band_path}: {e}")

def nbr(band1, band2):
    """
    Calculates the Normalized Burn Ratio (NBR) using the formula (NIR - SWIR) / (NIR + SWIR)
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        nbr = (band1 - band2) / (band1 + band2)
        nbr[np.isnan(nbr)] = 0  # Replace NaN values with 0
    return nbr

def dnbr(nbr1, nbr2):
    """
    Calculates dNBR by subtracting post-fire NBR from pre-fire NBR
    """
    return nbr1 - nbr2

def reproject_shp_gdal(infile, outfile, targetprj):
    """
    Reprojects a shapefile to match the projection of the satellite image (Sentinel-2)
    """
    try:
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(infile, 1)
        if dataSource is None:
            raise FileNotFoundError(f"Shapefile {infile} could not be opened.")

        layer = dataSource.GetLayer()
        sourceprj = layer.GetSpatialRef()
        transform = osr.CoordinateTransformation(sourceprj, targetprj)

        outDriver = ogr.GetDriverByName("Esri Shapefile")
        outDataSource = outDriver.CreateDataSource(outfile)
        outlayer = outDataSource.CreateLayer('', targetprj, ogr.wkbPolygon)
        outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))

        i = 0
        for feature in layer:
            transformed = feature.GetGeometryRef()
            transformed.Transform(transform)
            geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
            defn = outlayer.GetLayerDefn()
            feat = ogr.Feature(defn)
            feat.SetField('id', i)
            feat.SetGeometry(geom)
            outlayer.CreateFeature(feat)
            i += 1
            feat = None
    except Exception as e:
        raise RuntimeError(f"Error reprojecting shapefile {infile}: {e}")

def array2raster(array, geoTransform, projection, filename):
    """
    Converts a numpy array to a GeoTIFF format raster with given affine transformation and projection
    """
    try:
        pixels_x = array.shape[1]
        pixels_y = array.shape[0]

        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(
            filename,
            pixels_x,
            pixels_y,
            1,
            gdal.GDT_Float64)
        dataset.SetGeoTransform(geoTransform)
        dataset.SetProjection(projection)
        dataset.GetRasterBand(1).WriteArray(array)
        dataset.FlushCache()
        return dataset, dataset.GetRasterBand(1)
    except Exception as e:
        raise RuntimeError(f"Error creating raster file {filename}: {e}")

def clip_raster(filename, shp):
    """
    Clips a raster based on a shapefile using rasterio.
    Ensures that the CRS of the shapefile matches the raster.
    """
    try:
        inraster = rasterio.open(filename)
        
        # Ensure CRS match
        raster_crs = inraster.crs
        shapefile_crs = shp.crs
        
        if raster_crs != shapefile_crs:
            print(f"Warning: CRS mismatch! Reprojecting shapefile from {shapefile_crs} to {raster_crs}.")
            shp = shp.to_crs(raster_crs)  # Reproject shapefile to raster's CRS
            
        extent_geojson = [mapping(geom) for geom in shp.geometry]  # Ensure this is properly passed as a list of dicts
        
        clipped, crop_affine = mask(inraster, shapes=extent_geojson, nodata=np.nan, crop=True)
        
        clipped_meta = inraster.meta.copy()
        clipped_meta.update({
            "driver": "GTiff",
            "height": clipped.shape[1],
            "width": clipped.shape[2],
            "transform": crop_affine
        })
        cr_ext = rasterio.transform.array_bounds(clipped_meta['height'], clipped_meta['width'], clipped_meta['transform'])
        gt = crop_affine.to_gdal()
        return clipped, clipped_meta, cr_ext, gt
    except Exception as e:
        raise RuntimeError(f"Error clipping raster {filename}: {e}")


def reclassify(array):
    """
    Classifies the array into five severity levels based on dNBR values
    """
    reclass = np.zeros((array.shape[0], array.shape[1]))
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if np.isnan(array[i, j]):
                reclass[i, j] = np.nan
            elif array[i, j] < 0.1:
                reclass[i, j] = 1
            elif array[i, j] < 0.27:
                reclass[i, j] = 2
            elif array[i, j] < 0.44:
                reclass[i, j] = 3
            elif array[i, j] < 0.66:
                reclass[i, j] = 4
            else:
                reclass[i, j] = 5
    return reclass
