function emi_method (myFolder_tropomi, myFile_era5, inves_date, lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, pre_lev, sour_sou, allow_win_rot, prior_emi, EMG_method, GAU_method, CSF_method, plu_per)

[x_len, y_len, bin_lon_min, bin_lon_max, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat, inves_tropomi_data]  = data_extract (myFolder_tropomi, inves_date, lat_sou(1), lon_sou(1), inter_bin, tot_ran_y_abo, tot_ran_y_bel, tot_ran_x_upwind, tot_ran_x_dowwind, sour_sou); 

if ~isempty(inves_tropomi_data)

    meas_time_inves_mean = nanmean(inves_tropomi_data.meas_time);
    lat_inves = inves_tropomi_data.lat; lon_inves = inves_tropomi_data.lon; 
    no2_inves = inves_tropomi_data.no2;

    lat_bou1_inves = inves_tropomi_data.lat_bou1; lat_bou2_inves = inves_tropomi_data.lat_bou2;
    lat_bou3_inves = inves_tropomi_data.lat_bou3; lat_bou4_inves = inves_tropomi_data.lat_bou4;

    lon_bou1_inves = inves_tropomi_data.lon_bou1; lon_bou2_inves = inves_tropomi_data.lon_bou2;
    lon_bou3_inves = inves_tropomi_data.lon_bou3; lon_bou4_inves = inves_tropomi_data.lon_bou4;

    inves_date_str = num2str(inves_date); inves_year = inves_date_str(1:4);

    [era5_ws_ano, era5_wd_ano, era5_wd_ano_org, era5_ws_ano_uncer_1] = era5_specific_ano(myFile_era5, inves_year, meas_time_inves_mean, pre_lev, lat_sou(1), lon_sou(1));
    a_r = 180+era5_wd_ano-(allow_win_rot); 

    disp(['Wind speed (m/s): ', num2str(era5_ws_ano)])

    if ~isempty (era5_ws_ano)
        no2_bin = data_grid(lat_sou(1), lon_sou(1), inter_bin, tot_ran_y_abo, tot_ran_x_upwind, no2_inves, lat_inves, lon_inves, lat_bou1_inves, lat_bou2_inves, lat_bou3_inves, lat_bou4_inves, lon_bou1_inves, lon_bou2_inves, lon_bou3_inves, lon_bou4_inves, x_len, y_len, a_r, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat, plu_per);

        if contains(EMG_method ,'yes')
            [lifetime_est_emg, emission_est_emg, al, mu, x0, s, b] = emg_fitting(no2_bin, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano);
            disp('EMG method has been performed. Results are shown below:')
            disp(['Estimated lifetime (hr): ', num2str(lifetime_est_emg)])
            disp(['Estimated NO_X emission (g/s): ', num2str(emission_est_emg)])
            disp(['Estimated NO_X emission (Kt/year): ', num2str(emission_est_emg*10^(-9)*365*24*60*60)])
            disp(['Wind speed (m/s): ', num2str(era5_ws_ano)])
            disp(['Five fitted parameters: ', [' al: ', num2str(al)], [' mu: ', num2str(mu)], [' x0: ', num2str(x0)], [' s: ', num2str(s)], [' b: ', num2str(b)]])
            
            no2_bin_conv = no2_bin.* inter_bin*1000;
            lin_den = movmean(nansum(no2_bin_conv,1), 1);
            bkg_main = nanmean(lin_den (1:30/inter_bin), 'all');
            bkg_secon = nanmean(lin_den (1:20/inter_bin), 'all');

            bkg_uncer = (bkg_secon - bkg_main)/bkg_main;

            ws_uncer = (era5_ws_ano_uncer_1 - era5_ws_ano)/era5_ws_ano;

            total_uncer = sqrt((bkg_uncer^2)+ (ws_uncer^2));

            % disp(['Estimated NO_X emission uncertainity (Kt/year): ', num2str(emission_est_emg'*10^(-9)*365*24*60*60.*total_uncer)])

        end

        if contains(GAU_method ,'yes')
            

            [lifetime_est_gau, emission_est_gau] = gau_fitting(no2_bin, lat_sou, lon_sou, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano, prior_emi, a_r, tot_ran_x_dowwind, tot_ran_x_upwind);
            
            disp('Gau method has been performed. Results are shown below:')
            disp(['Estimated lifetime (hr): ', num2str(lifetime_est_gau')])
            disp(['Estimated NO_X emission (g/s): ', num2str(emission_est_gau')])
            disp(['Estimated NO_X emission (Kt/year): ', num2str(emission_est_gau'*10^(-9)*365*24*60*60)])
            
            no2_bin_conv = no2_bin.* inter_bin*1000;
            lin_den = movmean(nansum(no2_bin_conv,1), 1);
            bkg_main = nanmean(lin_den (1:30/inter_bin), 'all');
            bkg_secon = nanmean(lin_den (1:20/inter_bin), 'all');

            bkg_uncer = (bkg_secon - bkg_main)/bkg_main;

            ws_uncer = (era5_ws_ano_uncer_1 - era5_ws_ano)/era5_ws_ano;

            total_uncer = sqrt((bkg_uncer^2)+ (ws_uncer^2));

            % disp(['Estimated NO_X emission uncertainity (Kt/year): ', num2str(emission_est_gau'*10^(-9)*365*24*60*60.*total_uncer)])
        end

        if contains(CSF_method ,'yes')

            xData = (bin_lon_min + bin_lon_max)./2;
            xData_fit = xData([1:30/inter_bin, end-(30/inter_bin):end]);
            no2_bin_ave = nanmean(no2_bin, 1);
            no2_bin_ave_fit = no2_bin_ave([1:30/inter_bin, end-(30/inter_bin):end]);
            
            linear_fit = polyfit(xData_fit, no2_bin_ave_fit, 1);

            slope = linear_fit(1);
            intercept = linear_fit(2);

            linear_bkg_csf = slope*xData +intercept; 
            no2_bin_enh = (no2_bin - linear_bkg_csf) * 46.01 * inter_bin*1000;
            
            if ~exist('lifetime_est_gau', 'var')
                disp('Lifetime is not available from other methods. So user should provide the Lifetime for CSF method.');
                lifetime_est_gau = input('NOX lifetime (hr): ');
            end

            csf_emi_downwind = nan(1,3);

            for csf_dis_xx = 1:length(csf_emi_downwind)
                csf_emi = nansum(no2_bin_enh(:,(tot_ran_x_upwind/inter_bin)+csf_dis_xx), 'all') * era5_ws_ano;
                csf_emi_downwind(csf_dis_xx) = csf_emi/ exp(-(csf_dis_xx*(inter_bin)*1000)/(era5_ws_ano*lifetime_est_gau*3600));
            end

            disp('CSF has been performed. Results are shown below: ')
            disp(['Estimated NO_X emission (g/s): ', num2str(csf_emi_downwind)])
            disp(['Estimated mean NO_X emission (g/s): ', num2str(nanmean(csf_emi_downwind, 'all'))])
            disp(['Estimated NO_X emission (Kt/year): ', num2str(nanmean(csf_emi_downwind, 'all')*10^(-9)*365*24*60*60)])

            bkg_csf = nanmean(no2_bin(:,1:30/inter_bin), 'all');

            bkg_csf_uncer = nanmean(no2_bin(:,1:20/inter_bin), 'all');
            bkg_csf_uncer = (bkg_csf_uncer - bkg_csf)/bkg_csf;
            
            ws_uncer = (era5_ws_ano_uncer_1 - era5_ws_ano)/era5_ws_ano;

            total_uncer = sqrt((bkg_csf_uncer^2)+ (ws_uncer^2));

            disp(['Estimated NO_X emission uncertainity (Kt/year): ', num2str(nanmean(csf_emi_downwind, 'all')*10^(-9)*365*24*60*60 * total_uncer)])

        end

    elseif isempty (era5_ws_ano)
        disp('No ERA5 wind data is available for this power plant/investigation date. Program terminated.')
    end


elseif isempty(inves_tropomi_data)
    disp('No TROPOMI data is available for this power plant/investigation date. Program terminated.')
end

end