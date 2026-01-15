# Power Plant NO₂ Emission Estimation Using TROPOMI Data

This MATLAB project estimates **NO₂ emissions** from power plants using **TROPOMI satellite measurements** and **ERA5 meteorological data**.  
The workflow includes pre-processing satellite measurements, aligning emission plumes with wind direction, and estimating NO₂ emissions using one or more modeling methods.

---

## Overview

This script allows you to:  
- Extract and process TROPOMI NO₂ measurements over selected sources (e.g., power plants).  
- Align NO₂ plumes with wind direction from ERA5 data.  
- Identify and isolate emission plumes from nearby sources.  
- Estimate NOx **emission rates** and **atmospheric lifetime**.
- Choose among **three different methods** for emission estimation:  
  1. **Exponentially Modified Gaussian (EMG) method**  
  2. **Gaussian Plume Model (GPM) method**  
  3. **Cross-sectional Emission Flux (CSF) method**  

It supports **single or multiple emission sources** and allows customization of grid size, analysis domain, and plume detection settings.

---
## How to Run

- Modify the required input parameters in `Main_file.m`
- Run `Main_file.m`
---
## Required Input Files

| **Field** | **Description** | **Example / Notes** |
|-----------|----------------|------------------|
| `myFolder_tropomi` | Folder containing TROPOMI NO₂ data | `'D:\TROPOMI_NO2_india\data'` |
| `myFile_era5` | ERA5 meteorological data file (NetCDF) | `'D:\ERA5\ind_2022_7_9_900_1000.nc'` |
| `inves_date` | Date of investigation for the power plant | `20220224` (YYYYMMDD) |
| `lat_sou` | Latitude(s) of the emission source(s) | `16.499546` (can be array for multiple sources) |
| `lon_sou` | Longitude(s) of the emission source(s) | `75.834632` (can be array for multiple sources) |
| `prior_emi` | Prior emission estimate (g/s) | `repelem(100, length(lat_sou))` if unknown |
| `plu_per` | Percentile of NO₂ enhancement to define plume | `80` |
| `inter_bin` | Grid size in km | `5` → 5 km × 5 km grid |
| `tot_ran_y_abo` | Distance above source for emission calculation (km) | `40` |
| `tot_ran_y_bel` | Distance below source for emission calculation (km) | `40` |
| `tot_ran_x_upwind` | Upwind distance from source for background estimation (km) | `50` |
| `tot_ran_x_dowwind` | Downwind distance for emission calculation (km) | `80` |
| `pre_lev` | ERA5 pressure level(s) used for wind | `1000` (can be array, mean used) |
| `sour_sou` | Rectangular domain around source (degree) | `1` (no need to change) |
| `allow_win_rot` | Additional plume rotation angle (degrees) | `-10` (negative = clockwise) |
| `EMG_method` | Use EMG method? | `'yes'` or `'no'` |
| `GAU_method` | Use Gaussian plume method? | `'yes'` or `'no'` |
| `CSF_method` | Use cross-sectional fitting method? | `'yes'` or `'no'` |

---

