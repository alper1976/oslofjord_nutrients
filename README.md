# oslofjord_nutrients
# Climate-Driven Freshening and Nutrient Mitigation Drive Stoichiometric Decoupling and a Summer-to-Spring Regime Shift in the Oslofjord

This repository contains the data processing pipelines, statistical models, and visualization scripts required to reproduce the analyses and figures presented in the manuscript: 

> **Eiler, A., Puts, I., Stensrud, E., Werner, E., Ong, D., Langangen, Ø., Edvardsen, B., Santos, A. L. D., & Hessen, D. O. (in prep). Climate-Driven Freshening and Nutrient Mitigation Drive Stoichiometric Decoupling and a Summer-to-Spring Regime Shift in the Oslofjord.

## Overview
The provided R scripts synthesize nearly a century (1930–2025) of in situ monitoring data and modern satellite remote sensing to untangle the vertical decoupling of physical, chemical, and biological drivers of phytoplankton dynamics in the Oslofjord. 

Key analytical pipelines include:
* **Physical Forcing:** Calculating pycnocline deepening and stratification intensification using density gradients ($\Delta\rho/\Delta Z$).
* **Stoichiometry & GAMs:** Non-linear modeling of dissolved inorganic ($NO_3:PO_4$) and total (TN:TP) nutrient trajectories.
* **Multivariate Ordination:** Centered and scaled Principal Component Analysis (PCA) and Redundancy Analysis (RDA) demonstrating the vertical decoupling of the surface lens (salinity-driven) and the Deep Chlorophyll Maximum (stoichiometry-driven).
* **Phenodynamics:** Extracting spring bloom optimal temporal windows ("Goldilocks penalty") and modeling bimodal ecosystem regime shifts.

## Repository Structure
* `data/`: Contains both the raw inputs (from the Norwegian Environment Agency and GBIF) files required to run the scripts. 
* `scripts/`: Sequentially numbered R scripts covering data cleaning, statistical modeling, and plot generation.
* `outputs/`: Destination folder for the generated A4-formatted, publication-ready `.pdf` figures.

## Dependencies
The analyses were conducted in `R` (version 4.X.X). The following primary packages are required:

* **Data Wrangling & Spatial:** `tidyverse` (dplyr, tidyr, ggplot2), `lubridate`, `sf`
* **Statistical Modeling:** `mgcv` (GAMs), `vegan` (PCA/RDA)
* **Oceanography:** `oce` (Density calculations)
* **Visualization:** `factoextra`, `corrplot`, `patchwork`

To install the required packages, run:
```R
install.packages(c("tidyverse", "lubridate", "sf", "mgcv", "vegan", "oce", "factoextra", "corrplot", "patchwork"))
Reproduction Instructions
To reproduce the manuscript analyses:

Clone this repository to your local machine.


Run the scripts in the scripts/ folder sequentially, starting from XX

Figures 1 through 5, alongside Supplementary Figures S1–S8, will be automatically exported to the outputs/figures/ directory.

Data Availability
The foundational long-term monitoring data utilized in this study were aggregated from the Norwegian Environment Agency's Vannmiljø database and the GBIF Archive. Satellite-derived surface Chlorophyll-a data (MODIS) were acquired from the Copernicus Marine Environment Monitoring Service (CMEMS).

The compiled, spatially clustered datasets required to execute the scripts in this repository are permanently archived on Zenodo: [Insert Zenodo DOI Link Here].

License
The code in this repository is licensed under the MIT License - see the LICENSE file for details.