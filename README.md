# pyburn
Pyburn is a Python-based tool designed for calculating burn severity using Sentinel-2 imagery. With an emphasis on geospatial analysis, Pyburn leverages popular Python libraries like Rasterio, GeoPandas, and GDAL to provide robust functionality for raster processing and burn severity analysis.
graph TD
    A[Start] --> B[Read Pre-fire NIR and SWIR bands]
    B --> C[Read Post-fire NIR and SWIR bands]
    C --> D[Calculate Pre-fire NBR]
    D --> E[Calculate Post-fire NBR]
    E --> F[Calculate dNBR - Post-fire NBR minus Pre-fire NBR]
    F --> G[Visualize Pre-fire NBR]
    G --> H[Visualize Post-fire NBR]
    H --> I[Visualize dNBR]
    I --> J[Reproject Shapefile to Match Sentinel-2 Projection]
    J --> K[Clip dNBR Using Fire Boundary Shapefile]
    K --> L[Save dNBR to GeoTIFF]
    L --> M[Plot Burn Severity Map Using Reclassified dNBR]
    M --> N[Calculate Burnt Area by Severity Class]
    N --> O[End]
