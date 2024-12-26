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


# Paths to Sentinel-2 images (pre-fire and post-fire)
path_prefire_nir = r"C:\Users\Dell\Downloads\Pre fire satellite images sentinel 2\S2B_MSIL2A_20210728T094029_N0500_R036_T33SXD_20230220T230541.SAFE\GRANULE\L2A_T33SXD_A022942_20210728T094408\IMG_DATA\R20m\T33SXD_20210728T094029_B8A_20m.jp2"
path_prefire_swir = r"C:\Users\Dell\Downloads\Pre fire satellite images sentinel 2\S2B_MSIL2A_20210728T094029_N0500_R036_T33SXD_20230220T230541.SAFE\GRANULE\L2A_T33SXD_A022942_20210728T094408\IMG_DATA\R20m\T33SXD_20210728T094029_B12_20m.jp2"
path_postfire_nir = r"C:\Users\Dell\Downloads\Post fire satellite images sentinel 2\S2B_MSIL2A_20210817T094029_N0500_R036_T33SXD_20230216T200808.SAFE\GRANULE\L2A_T33SXD_A023228_20210817T094030\IMG_DATA\R20m\T33SXD_20210817T094029_B8A_20m.jp2"
path_postfire_swir = r"C:\Users\Dell\Downloads\Post fire satellite images sentinel 2\S2B_MSIL2A_20210817T094029_N0500_R036_T33SXD_20230216T200808.SAFE\GRANULE\L2A_T33SXD_A023228_20210817T094030\IMG_DATA\R20m\T33SXD_20210817T094029_B12_20m.jp2"
# Path to shapefile
infile_shp = r"C:\Users\Dell\Downloads\crotene zip\crotene.shp"

# Path for the reprojected shapefile
outfile_shp = r"C:\Users\Dell\Downloads\reprojected_shapefile.shp"  # save the reprojected_shapefile in your preferred file path 

# Output filenames
filename = r"C:\Users\Dell\Downloads\output_dnbr.tiff"   # save the dnbr in your preferred file path 
filename2 = r"C:\Users\Dell\Downloads\output_clipped_dnbr.tiff"   # save the clipped_dnbr in your preferred file path
fname = r"C:\Users\Dell\Downloads\output_burn_severity_map.png"    # save the burn_severity_map in your preferred file path


## Read pre-fire NIR and SWIR images
pre_fire_nir, crs, geoTransform, targetprj = read_band_image(path_prefire_nir)
pre_fire_swir, crs, geoTransform, targetprj = read_band_image(path_prefire_swir)

# Read post-fire NIR and SWIR images
post_fire_nir, crs, geoTransform, targetprj = read_band_image(path_postfire_nir)
post_fire_swir, crs, geoTransform, targetprj = read_band_image(path_postfire_swir)

# Calculate pre-fire NBR
pre_fire_nbr = nbr(pre_fire_nir.astype(int), pre_fire_swir.astype(int))

# Calculate post-fire NBR
post_fire_nbr = nbr(post_fire_nir.astype(int), post_fire_swir.astype(int))

# Calculate dNBR (pre-fire NBR - post-fire NBR)
DNBR = dnbr(pre_fire_nbr, post_fire_nbr)


# VISUALIZE THE IMAGES
# ---------------------------------------------------
# 1. Visualize Pre-fire NBR
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'xticks': [], 'yticks': []})
cax = ax.imshow(pre_fire_nbr, cmap='RdYlGn', vmin=-0.5, vmax=1.0)
plt.title('Pre-fire NBR')
cbar = fig.colorbar(cax, ax=ax, fraction=0.035, pad=0.04)
plt.show()
plt.savefig('pre_fire_nbr.png', bbox_inches="tight")

# 2. Visualize Post-fire NBR
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'xticks': [], 'yticks': []})
cax = ax.imshow(post_fire_nbr, cmap='RdYlGn', vmin=-0.5, vmax=1.0)
plt.title('Post-fire NBR')
cbar = fig.colorbar(cax, ax=ax, fraction=0.035, pad=0.04)
plt.show()
plt.savefig('post_fire_nbr.png', bbox_inches="tight")

# 3. Visualize dNBR (Post-fire NBR - Pre-fire NBR)
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'xticks': [], 'yticks': []})
cax = ax.imshow(DNBR, cmap='RdYlGn', vmin=-0.5, vmax=1.0)
plt.title('dNBR (Post-fire - Pre-fire)')
cbar = fig.colorbar(cax, ax=ax, fraction=0.035, pad=0.04)
plt.show()
plt.savefig('dnbr.png', bbox_inches="tight")
# ---------------------------------------------------


# Reproject the shapefile to match the Sentinel-2 image projection
reproject_shp_gdal(infile_shp, outfile_shp, targetprj)

# Read the reprojected shapefile
fire_boundary = gpd.read_file(outfile_shp)

# Save dNBR to a GeoTIFF
dnbr_tif, dnbr_tifBand = array2raster(DNBR, geoTransform, crs, filename)

# Clip the dNBR raster using the reprojected shapefile
(clipped_dnbr, clipped_dnbr_meta, cr_extent, gt) = clip_raster(filename, fire_boundary)

# Save the clipped dNBR to a new GeoTIFF
clipped_ds, clipped_ds_rasterband = array2raster(clipped_dnbr[0], gt, crs, filename2)

# Plot the burn severity map (which you already have)
# Custom burn severity colors
from matplotlib.colors import ListedColormap

# Define your custom colormap using the specified colors: green, yellow, orange, red, purple
cmap = ListedColormap(['green', 'yellow', 'orange', 'red', 'purple'])

# Normalize to the severity levels 1 to 5 (based on reclassification)
norm = Normalize(vmin=1, vmax=5)  # Values corresponding to reclassification (1 to 5)

# Plot the burn severity map based on reclassified values (not raw dNBR)
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'xticks': [], 'yticks': []})

# Use the reclassified array for plotting
reclassified_array = reclassify(clipped_ds_rasterband.ReadAsArray())

# Display the reclassified burn severity map
cax = ax.imshow(reclassified_array, cmap=cmap, norm=norm)

# Title and colorbar
plt.title('Burn Severity Map')
cbar = fig.colorbar(cax, ax=ax, fraction=0.035, pad=0.04, ticks=[1, 2, 3, 4, 5])
cbar.ax.set_yticklabels(['Unburned', 'Low Severity', 'Moderate-low Severity', 'Moderate-high Severity', 'High Severity'])

# Show and save the figure
plt.show()
plt.savefig(fname, bbox_inches="tight")


# Calculate burnt area based on reclassified dNBR
reclass = reclassify(clipped_ds_rasterband.ReadAsArray())
severity_labels = ['Unburned hectares', 'Low severity hectares', 'Moderate-low severity hectares', 'Moderate-high severity hectares', 'High severity']
for i in range(1, 6):
    x = reclass[reclass == i]
    area = x.size * 0.04  # Pixel size (20m x 20m)
    print(f"{severity_labels[i-1]}: {area:.2f} hectares")