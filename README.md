# MEARDN Sandfish Density Analysis

This repository contains data and analysis code for the manuscript:
“Conservation and Commerce: Community-Based Management of Sandfish in Papua New Guinea”.

## Contents
- MEARDN_Sandfish_UVC_Data.csv
Underwater visual census (UVC) survey data on sandfish (Holothuria scabra) collected by community monitors between 2017–2019.
Fields include: site code, date, observer ID, tenure status (LMMA or open), and sandfish count.

- TenureAreaCombinationPlots_V2.m
Annotated MATLAB script used to process the data and generate violin plots of sandfish densities by site, year, and tenure. This script reproduces Figure 1 in the main manuscript.

## Usage
1. Requirements:
- MATLAB R2021a or later (developed in MATLAB R2024b).
- No external dependencies required.
2. Running the analysis:
- Place both files in the same directory.
- Open TenureAreaCombinationPlots_V2.m in MATLAB.
- Adjust file names/paths at the top of the script if needed.
- Run the script to produce summary statistics and generate violin plots (saved as PDF and PNG).
3. Output:
- Publication-quality figures comparing sandfish densities inside and outside LMMAs across three survey periods.
- A CSV file with summary statistics (mean densities and bootstrap confidence intervals for each group).

## Data Sensitivity
- Site coordinates have been anonymised to protect the locations of locally managed marine areas (LMMAs).

## Citation

If you use these data or scripts, please cite:
Hamilton, R. et al. “Conservation and Commerce: Community-Based Management of Sandfish in Papua New Guinea” (TBD).

## Contact

For questions, please contact michael.bode@qut.edu.au
