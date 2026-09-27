# Urbanization destabilizes rhizosphere multitrophic networks with cascading risks to urban One Health

This repository contains the reproducible R scripts and analysis pipelines for the manuscript:
**"Urbanization destabilizes rhizosphere multitrophic networks with cascading risks to urban One Health"**

## Overview of R Scripts
All analysis workflows correspond directly to the figures presented in the manuscript:
- `01_Fig2_Biogeochemistry_Suite.R`: Biogeochemical profiling, heavy metal risk modeling (PLI, RI, HI), and ecoenzymatic vector stoichiometry.
- `02_Fig3_Diversity_DDR.R`: Tri-kingdom alpha/beta diversity, taxonomic replacement vs. richness difference partitioning, and distance-decay relationships (DDR).
- `03_Fig4_Niche_Phylogeny.R`: Multitrophic niche breadths, Blomberg's K phylogenetic conservation, and community adaptability.
- `04_Fig5_Assembly_Processes.R`: Phylogenetic null models (βNTI, RCbray), Sloan neutral community models (NCM), and environmental sensitivity classification.
- `05_Fig6_Multifunctionality.R`: Cross-trophic co-occurrence networks, in silico topological robustness simulations, and EMF decoupling.
- `06_Fig7_Fig8_SoilHealth_OneHealth.R`: CASH soil health scoring, Safe Operating Space (SOS) modeling, trophic optimum (1.64), and Structural Equation Modeling (SEM).

## System Requirements & Dependencies
- **Software**: R (version >= 4.3.0) / RStudio
- **Operating Systems**: Windows 10/11, macOS, or Linux
- **Required R packages**:
  - Community ecology: `vegan`, `picante`, `linkET`, `iNEXT`
  - Modeling & Statistics: `lme4`, `lavaan`, `smatr`
  - Networks & Visualization: `igraph`, `ggplot2`, `patchwork`, `readxl`

## Instructions for Use
1. Clone or download this repository to your local machine.
2. Open R / RStudio and set the working directory to the folder containing these scripts:
   ```R
   setwd("path/to/Rhizosphere-Urban-OneHealth")
