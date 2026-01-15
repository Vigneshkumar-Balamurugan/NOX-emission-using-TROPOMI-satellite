
function [era5_ws_ano, era5_wd_ano, era5_wd_ano_org, era5_ws_ano_uncer_1] =  era5_specific_ano(myFile_era5, inves_year, meas_time_inves_mean, pre_lev, lat_sou, lon_sou)


myFile_era5 = fullfile(myFolder_era5, era5_file_name);

era5_file_info = ncinfo(myFile_era5);

era5_u = ncread(myFile_era5,'u'); era5_v = ncread(myFile_era5,'v'); 
era5_lat = ncread(myFile_era5,'latitude'); era5_lon = ncread(myFile_era5,'longitude'); 

if ~isempty(find(strcmp({era5_file_info.Variables.Name}, 'level')))
    era5_level = ncread(myFile_era5,'level');
else
    era5_level = ncread(myFile_era5,'pressure_level');
end

if ~isempty(find(strcmp({era5_file_info.Variables.Name}, 'time')))
    era5_time = ncread(myFile_era5,'time'); 
    time_mani = datetime(1970, 1, 1, 0, 0, 0, 'TimeZone','UTC') - datetime(1900, 1, 1, 0, 0, 0, 'TimeZone','UTC');
    era5_time_cha = (era5_time - hours(time_mani))*60*60; % epoch time (seconds since 1970)
else
    era5_time = ncread(myFile_era5,'valid_time'); % epoch time (seconds since 1970)
    era5_time_cha = era5_time;
end



    
[~, era5_lat_near_ano] = (min(abs(lat_sou-era5_lat))); 
[~, era5_lon_near_ano] = (min(abs(lon_sou-era5_lon)));

[~, era5_time_near_ano] = (min(abs(meas_time_inves_mean-era5_time_cha)));

if abs(min(abs(meas_time_inves_mean-era5_time_cha))) <= 10800
    era5_u_ano = (era5_u(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano-1:era5_time_near_ano+1)); 
    era5_v_ano = (era5_v(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano-1:era5_time_near_ano+1));
    
    era5_u_ano_inter = interp1(double(era5_time_cha(era5_time_near_ano-1:era5_time_near_ano+1)), squeeze(nanmean(era5_u_ano, 3)), meas_time_inves_mean, 'linear');
    era5_v_ano_inter = interp1(double(era5_time_cha(era5_time_near_ano-1:era5_time_near_ano+1)), squeeze(nanmean(era5_v_ano, 3)), meas_time_inves_mean, 'linear');

    era5_u_ano_aprox = nanmean((era5_u(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano)), 'all'); 
    era5_v_ano_aprox = nanmean((era5_v(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, pre_lev)),era5_time_near_ano)), 'all');

    era5_ws_ano = ((era5_u_ano_inter.^2) + (era5_v_ano_inter.^2)).^(1/2); 
    era5_wd_ano_org = 180 + ((180/ pi) * atan2 (era5_u_ano_aprox, era5_v_ano_aprox));

    if era5_wd_ano_org >= 0 && era5_wd_ano_org <= 90
        era5_wd_ano = 270 - era5_wd_ano_org ;
    elseif era5_wd_ano_org >= 90 && era5_wd_ano_org <= 180
        era5_wd_ano = 180 - era5_wd_ano_org + 90;
    elseif era5_wd_ano_org >= 180 && era5_wd_ano_org <= 270
    era5_wd_ano = 270 - era5_wd_ano_org;
    elseif era5_wd_ano_org >= 270 && era5_wd_ano_org <= 360
        era5_wd_ano = 360 - era5_wd_ano_org + 270;
    end
    
    era5_u_ano_uncer_1 = (era5_u(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, 975)),era5_time_near_ano-1:era5_time_near_ano+1)); 
    era5_v_ano_uncer_1 = (era5_v(era5_lon_near_ano,era5_lat_near_ano,find(ismember(era5_level, 975)),era5_time_near_ano-1:era5_time_near_ano+1));

    era5_u_ano_inter_uncer_1 = interp1(double(era5_time_cha(era5_time_near_ano-1:era5_time_near_ano+1)), squeeze(nanmean(era5_u_ano_uncer_1, 3)), meas_time_inves_mean, 'linear');
    era5_v_ano_inter_uncer_1 = interp1(double(era5_time_cha(era5_time_near_ano-1:era5_time_near_ano+1)), squeeze(nanmean(era5_v_ano_uncer_1, 3)), meas_time_inves_mean, 'linear');

    era5_ws_ano_uncer_1 = ((era5_u_ano_inter_uncer_1.^2) + (era5_v_ano_inter_uncer_1.^2)).^(1/2); 
else 
    era5_ws_ano = [];
    era5_wd_ano = [];
    era5_wd_ano_org = [];

    era5_ws_ano_uncer_1 = [];
end

end
