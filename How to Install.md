# PyBurn: Geospatial Analysis for Burn Severity

PyBurn is a Python-based geospatial library for calculating burn severity using Sentinel-2 satellite images. The library calculates Normalized Burn Ratio (NBR), differential NBR (dNBR), visualizes pre-fire and post-fire conditions, and classifies burn severity into various levels. This project also includes an analysis of burnt areas.

### Features
- Reads and processes Sentinel-2 satellite images (NIR and SWIR bands).
- Calculates Pre-fire and Post-fire NBR and dNBR.
- Visualizes NBR and dNBR images.
- Reprojects shapefiles to match Sentinel-2 satellite image projection.
- Clips rasters based on fire boundary shapefiles.
- Classifies burn severity into 5 levels.
- Calculates the burnt area for each severity class in hectares.

---

## Prerequisites

- **Python** (version 3.x)
- **Anaconda** (for creating and managing environments)

Due to the large size of the Sentinel-2 satellite images, they are **not hosted on GitHub**. You can download the images and shapefiles from **Google Drive** and extract them before setting up the environment.

---

## How to Use PyBurn

### 1. Clone or Download the Repository

You can either clone the repository using Git or directly download the folder.

#### Option 1: Clone the Repository
#### Option 2: Download the ZIP
You can download the repository as a ZIP file from GitHub. Simply extract the contents to your local machine.
https://drive.google.com/file/d/1x6e2aCE64jAOcS_zZAVkXIFd1jUvSaIG/view?usp=drive_link

### 2. Download Satellite Images and Shapefiles
Due to the large size of the satellite images, they are hosted on Google Drive. You can download the images and shapefiles from the link below:
https://drive.google.com/file/d/1K-8Rgf0WHt_eyOMtW87YU5sEY83FowDe/view?usp=drive_link - satellite images link
https://drive.google.com/file/d/16qLa7CjLTo5YY_6YaBo_b3VtSTSakJ6Z/view?usp=drive_link - shapefile link

The zip folder contains:

Pre-fire NIR (Band 8A) and SWIR (Band 12) Sentinel-2 images.
Post-fire NIR (Band 8A) and SWIR (Band 12) Sentinel-2 images.
Shapefile folder for fire boundary.
After downloading, extract the files to a local folder.

### 3. Install Dependencies
To set up the environment and install the necessary dependencies, follow the steps below:

#### Step 1: Create the Conda Environment
Ensure that you have Anaconda installed. Then, in the Anaconda Prompt, navigate to the folder where environment.yml is located and run:
#### conda env create -f environment.yml

#### Step 2: Activate the Conda Environment
Activate the newly created environment:
#### conda activate gsp_pyburn_env

#### Step 3: Install the PyBurn Library
Now install the PyBurn library:
#### pip install .

### 4. Configure Paths in pyburn.py
After installation, you need to update the paths to the satellite images and shapefiles in the jupyter lab while testing with the data.

Pre-fire NIR and SWIR paths
Update the paths to the pre-fire NIR and SWIR images.

Post-fire NIR and SWIR paths
Update the paths to the post-fire NIR and SWIR images.

Shapefile Path
Update the path to the shapefile folder for the fire boundary.

Output Paths
Set the paths for the output files where you want to save the processed files.

Example Configuration
In jupyter lab, update the following variables:

# Paths to Sentinel-2 images (pre-fire and post-fire)

#### path_prefire_nir = r"C:\path\to\your\pre_fire_nir_image.jp2"

#### path_prefire_swir = r"C:\path\to\your\pre_fire_swir_image.jp2"

#### path_postfire_nir = r"C:\path\to\your\post_fire_nir_image.jp2"

#### path_postfire_swir = r"C:\path\to\your\post_fire_swir_image.jp2"

# Path to shapefile

#### infile_shp = r"C:\path\to\your\shapefile\crotene.shp"

# Path for the reprojected shapefile

#### outfile_shp = r"C:\path\to\save\reprojected_shapefile.shp"

# Output filenames

#### filename = r"C:\path\to\save\output_dnbr.tiff"

#### filename2 = r"C:\path\to\save\output_clipped_dnbr.tiff"

#### fname = r"C:\path\to\save\burn_severity_map.png"

### 5. Running the Script
You can also test the library functions of pyburn.py by running the test_pyburn.py script in your Anaconda Prompt by navigating to the directory containing the script and executing:
#### python -m unittest test_pyburn.py -v

As we can see all the below functions will run correctly:

test_array2raster (test_pyburn.TestPyburn) ... ok

test_clip_raster (test_pyburn.TestPyburn) ... ok

test_dnbr (test_pyburn.TestPyburn)
Test the dNBR calculation. ... ok

test_nbr (test_pyburn.TestPyburn)
Test the NBR calculation. ... ok

test_read_band_image (test_pyburn.TestPyburn) ... ok

test_reclassify (test_pyburn.TestPyburn) ... ok

test_reproject_shp_gdal (test_pyburn.TestPyburn)
Test the reproject_shp_gdal function. ... ok

----------------------------------------------------------------------
Ran 7 tests in 0.458s

OK
