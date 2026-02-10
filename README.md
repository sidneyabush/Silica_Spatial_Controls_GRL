# Thinking outside the rocks: Subsurface water storage, topography, and land cover are key modulators of large‐scale riverine dissolved silicon dynamics

**Bush, S. A.**, Johnson, K., Jankowski, K. J., Carey, J. C., Sethna, L. R., Lyon, N. J., & Sullivan, P. L. (2026). Thinking outside the rocks: Subsurface water storage, topography, and land cover are key modulators of large‐scale riverine dissolved silicon dynamics. Geophysical Research Letters, 53(2), e2025GL118853.  https://doi.org/10.1029/2025GL118853


## Data availability 

**Primary source data.** All river chemistry and ancillary data used here originate from the following USGS data release (please cite if you reuse):

> Jankowski, K. J., Johnson, K., Carey, J. C., Lyon, N., Julian, P., Bush, S., Sethna, L. R., Chen, A., Wymore, A. S., Kortelainen, P., Laudon, H., Poste, A., McKnight, D. M., McDowell, W. H., Shogren, A. J., Heindel, R. C., Raike, A., Jones, J. B., & Sullivan, P. L. (2025). **Global Aggregation of Stream Silica (GlASS) (Version 2.0, July 2025)** [Data release]. U.S. Geological Survey. https://doi.org/10.5066/P138M8AR

**Analysis-ready inputs.** To facilitate exact reproduction of the figures and results, the harmonized/partitioned inputs generated from the USGS release (e.g., driver tables and split definitions) are archived on Zenodo: 
* **Data: [![DOI](https://zenodo.org/badge/DOI/105281/zenodo.16884281.svg)](https://doi.org/10.5281/zenodo.16884281)**
* **Software: [![DOI](https://zenodo.org/badge/1038792577.svg)](https://doi.org/10.5281/zenodo.16884258)**  

> **Notes:**  
> - The Step 1 scripts document the transformations from the USGS release to analysis-ready tables; see each script’s header for expected filenames and paths.


## Repository structure

- `Step1_Harmonization/` — scripts to convert units, assemble WRTDS–Kalman + discharge tables for Catalina Jemez sites, and build the harmonized drivers + data partitions for all sites.
- `Step2_RF_Model_SHAP/` — scripts to train RF models (FNConc, FNYield), generate predictions/diagnostics, and compute SHAP values.
- `Step3_Create_Publication_Figures/` — scripts to build all manuscript and SI figures.


## Software & R packages

- **R ≥ 4.3.x** on macOS  
- Core packages used across scripts include:  
  `data.table`, `dplyr`, `tidyr`, `stringr`, `lubridate`,  
  `ggplot2`, `cowplot`, `patchwork`, `sf`, `ggspatial`, `ggrepel`, `maps`,  
  `randomForest`, `fastshap`, `iml`, `pdp`,  
  `foreach`, `doParallel`, and `librarian`.

