# PyBurn Script Flowchart

## Process Overview

The script follows a sequence of steps to process Sentinel-2 satellite images (pre-fire and post-fire) for calculating burn severity based on the Normalized Burn Ratio (NBR) and dNBR (differential NBR). It visualizes the results and calculates the burnt area in hectares.

---

### Flowchart

```mermaid
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
