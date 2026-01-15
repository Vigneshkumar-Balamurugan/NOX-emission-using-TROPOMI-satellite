function [lifetime_gau, emission_gau] = gau_fitting(no2_bin, lat_sou, lon_sou, inter_bin, bin_lon_min, bin_lon_max, era5_ws_ano, prior_emi, a_r, tot_ran_x_dowwind, tot_ran_x_upwind)

% fitting the unknown parameters, and estimating NOX emission and lifetime
no2_bin_conv = no2_bin.* inter_bin*1000;
lin_den = movmean(nansum(no2_bin_conv,1), 1);

bkg = nanmean(lin_den (1:30/inter_bin), 'all'); % first 30km will be used for calculating Background line density
yData = lin_den - bkg;

xData = (bin_lon_min + bin_lon_max)./2;


% lower and upper bound of 2 unknown fitted parameters
lb = [zeros(1, length(prior_emi)) 0]; 
ub = [repelem(10000, length(prior_emi)) 24]; 

lb = [zeros(1, length(prior_emi)) zeros(1, length(prior_emi))]; 
ub = [repelem(inf, length(prior_emi)) repelem(24, length(prior_emi))]; 


% Initial parameter guess
initialParams = [prior_emi, repelem(4, length(prior_emi))];

% Options for fmincon
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

% Perform the optimization
params = fmincon(@(params) gaufittingObjective(params, yData, lat_sou, lon_sou, era5_ws_ano, a_r-180, tot_ran_x_dowwind, tot_ran_x_upwind, inter_bin, prior_emi), initialParams, [], [], [], [], lb, ub, [], options);

% fitted parameters
emission_gau = [];
lifetime_gau = [];
for parm_len = 1:length(prior_emi)
    emission_gau = [emission_gau; params(parm_len)];
    lifetime_gau = [lifetime_gau; params(parm_len+length(prior_emi))];
end


[gau_model_bin, along_wind, no2_vc_summed] =  gau_model(emission_gau, lifetime_gau, lat_sou, lon_sou, era5_ws_ano, a_r-180);


no2_ld = nansum(no2_vc_summed,1);
no2_ld_req = no2_ld(((length(along_wind)-1)/2)+2:find(along_wind == tot_ran_x_dowwind));

no2_ld_req_mean =  nanmean(reshape(no2_ld_req, [inter_bin/gau_model_bin, tot_ran_x_dowwind/inter_bin]),1);

model = zeros(1,(tot_ran_x_upwind+tot_ran_x_dowwind)/inter_bin);

model(tot_ran_x_upwind/inter_bin+1:(tot_ran_x_upwind/inter_bin)+tot_ran_x_dowwind/inter_bin) = no2_ld_req_mean;

yData_2 = yData;

figure     
plot (xData, yData_2, 'color','b','linewidth',3,'LineStyle','-','DisplayName', 'Observed')
hold on
plot (xData, model, 'color','r','linewidth',3,'LineStyle','-','DisplayName', 'Fitted')
grid minor
legend()
ylim([0 max(yData)+3])
xlabel('distance (km)')
ylabel('Line density (mole/m)')
set(gca, 'YDir', 'normal','FontSize', 14,'fontweight','bold','FontName', 'Times New Roman')

figure 
imagesc(no2_vc_summed)

end