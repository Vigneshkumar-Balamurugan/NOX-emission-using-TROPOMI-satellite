% This script estimates point-source NOx emissions using TROPOMI satellite observations using Gaussian Plume model
clear; clc; close all;
%% Required Fields

myFolder_tropomi = 'D:\TROPOMI_NO2_india\data'; % TROPOMI NO2 folder
myFile_era5 = 'D:\ERA5\ind_2022_7_9_900_1000.nc';  % ERA5 file

inves_date = 20220224; % Investigation date of the power plant (in YYYYMMDD)

lat_sou = 16.499546; 
lon_sou = 75.834632;
% lat_sou = [20.90056, 21.1427]; % E.g., for for multiple sources 
% lon_sou = [85.1961, 85.0839];

prior_emi = repelem(100, length(lat_sou)); % in Prior emission in g/s (use 100 if unknown)

plu_per = 80; % % Percentile of enhancement considered as a plume

inter_bin = 5;            % Grid size (1 = 1 km × 1 km)
tot_ran_y_abo = 40;       % Distance (km) above the emission source used in emission calculation
tot_ran_y_bel = 40;       % Distance (km) below the emission source used in emission calculation
tot_ran_x_upwind = 50;    % Distance (km) in the upwind direction of the emission source used in emission calculation
tot_ran_x_dowwind = 80;  % Distance (km) in the downwind direction of the emission source used in emission calculation

pre_lev = 1000;  % ERA5 pressure level used for wind information in this analysis and emission estimation
                 % (can be multiple levels, e.g., [1000, 975, 950]; the mean value will be used)

sour_sou = 1;    % Rectangular domain around the source considered in this analysis (1 = 1 degree). No need to change.

allow_win_rot = -10; % Angle to rotate the plume in addition to ERA5 wind info
                     % Negative values rotate clockwise; positive values rotate counterclockwise
                     

EMG_method = 'yes'; GAU_method = 'yes'; CSF_method = 'no'; % Select "yes" or "no" to indicate whether to use the specific method.
%% Emission Estimation: Exponentially Modified Gaussian method, Gaussian Plume Model and/or Cross-sectional Emission Flux Method

emi_method (myFolder_tropomi, myFile_era5, inves_date, lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, pre_lev, sour_sou, allow_win_rot, prior_emi, EMG_method, GAU_method, CSF_method, plu_per);


