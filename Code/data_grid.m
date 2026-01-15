function no2_bin = data_grid(lat_sou, lon_sou, inter_bin, tot_ran_y_abo, tot_ran_x_upwind, no2_inves, lat_inves, lon_inves, lat_bou1_inves, lat_bou2_inves, lat_bou3_inves, lat_bou4_inves, lon_bou1_inves, lon_bou2_inves, lon_bou3_inves, lon_bou4_inves, x_len, y_len, a_r, bin_lat_min_mat, bin_lat_max_mat, bin_lon_min_mat, bin_lon_max_mat, bin_lon_mat, bin_lat_mat, plu_per)

% Data manipulation (i.e, TROPOMI NO2 rotation and gridding)

cen_xdis = nan(1,length(lat_inves)); cen_ydis = nan(1,length(lon_inves));

bou1_xdis = nan(1,length(lat_inves)); bou2_xdis = nan(1,length(lat_inves));
bou3_xdis = nan(1,length(lat_inves)); bou4_xdis = nan(1,length(lat_inves));
        
bou1_ydis = nan(1,length(lat_inves)); bou2_ydis = nan(1,length(lat_inves));
bou3_ydis = nan(1,length(lat_inves)); bou4_ydis = nan(1,length(lat_inves));
    
for no_ele = 1:length(lon_inves)
    
    [cen_xdis(no_ele), cen_ydis(no_ele)] = diskm (a_r, lat_sou, lon_sou, lat_inves(no_ele), lon_inves(no_ele));
    
    [bou1_xdis(no_ele), bou1_ydis(no_ele)] = diskm (a_r, lat_sou, lon_sou, lat_bou1_inves(no_ele), lon_bou1_inves(no_ele));
    [bou2_xdis(no_ele), bou2_ydis(no_ele)] = diskm (a_r, lat_sou, lon_sou, lat_bou2_inves(no_ele), lon_bou2_inves(no_ele));
    [bou3_xdis(no_ele), bou3_ydis(no_ele)] = diskm (a_r, lat_sou, lon_sou, lat_bou3_inves(no_ele), lon_bou3_inves(no_ele));
    [bou4_xdis(no_ele), bou4_ydis(no_ele)] = diskm (a_r, lat_sou, lon_sou, lat_bou4_inves(no_ele), lon_bou4_inves(no_ele));       
end

bou_xdis = [bou1_xdis; bou2_xdis; bou3_xdis; bou4_xdis];
bou_ydis = [bou1_ydis; bou2_ydis; bou3_ydis; bou4_ydis];

no2_bin = nan(y_len,x_len);

for y_len_xx = 1:y_len
    for x_len_xx = 1:x_len
        bin_loop_poly = polyshape([bin_lon_min_mat(y_len_xx,x_len_xx),bin_lat_min_mat(y_len_xx,x_len_xx);bin_lon_max_mat(y_len_xx,x_len_xx),bin_lat_min_mat(y_len_xx,x_len_xx);bin_lon_max_mat(y_len_xx,x_len_xx),bin_lat_max_mat(y_len_xx,x_len_xx);bin_lon_min_mat(y_len_xx,x_len_xx),bin_lat_max_mat(y_len_xx,x_len_xx)]); 
        
        [~, near_pix_idx] = find ((cen_xdis <= (bin_lon_max_mat(y_len_xx,x_len_xx)+3) & cen_xdis >= (bin_lon_min_mat(y_len_xx,x_len_xx)-3))...
                & (cen_ydis <= (bin_lat_max_mat(y_len_xx,x_len_xx)+2) & cen_ydis >= (bin_lat_min_mat(y_len_xx,x_len_xx)-2)));
        
        no2_loop_val = [];
        for pix_idx = 1:length(near_pix_idx)
            bou_loop_poly = polyshape([bou_xdis(1,near_pix_idx(pix_idx)),bou_ydis(1,near_pix_idx(pix_idx));bou_xdis(2,near_pix_idx(pix_idx)),bou_ydis(2,near_pix_idx(pix_idx));bou_xdis(3,near_pix_idx(pix_idx)),bou_ydis(3,near_pix_idx(pix_idx));bou_xdis(4,near_pix_idx(pix_idx)),bou_ydis(4,near_pix_idx(pix_idx))]);
            int_loop = intersect(bin_loop_poly, bou_loop_poly);
            
            if ~isempty(int_loop.Vertices)
                int_loop_area = polyarea(int_loop.Vertices(:,1),int_loop.Vertices(:,2));
                wei_area_loop = int_loop_area/(inter_bin*inter_bin);
                no2_loop_val = [no2_loop_val;no2_inves(near_pix_idx(pix_idx))*wei_area_loop,wei_area_loop];
            end
        end
        
        if ~isempty(no2_loop_val)
            no2_bin (y_len_xx, x_len_xx) = sum(no2_loop_val(:,1),'all') / sum(no2_loop_val(:,2),'all');
        end
    
    end
end


figure
scatter(lon_inves, lat_inves, [], no2_inves, 'o', 'filled')
grid minor
xlim([min(lon_inves(:)), max(lon_inves(:))])
ylim([min(lat_inves(:)), max(lat_inves(:))])
xlabel('longitude')
ylabel('latitude')
title('TROPOMI data')
set(gca,'FontSize', 15,'fontweight','bold','FontName', 'Times New Roman')

figure
ima_bin = imagesc(bin_lon_mat(1,:), bin_lat_mat(:,1), no2_bin);
set(ima_bin, 'AlphaData', ~isnan(no2_bin))
xlabel('distance (km)')
ylabel('distance (km)')
title('TROPOMI data gridded')
set(gca, 'YDir', 'normal','FontSize', 12,'fontweight','bold','FontName', 'Times New Roman')
pbaspect([4 1 1])
axis tight;

bkg_value = prctile(no2_bin(:), 25);

mask = isnan(no2_bin);
no2_bin = regionfill(no2_bin, mask);
no2_bin = fillmissing(no2_bin, 'constant', bkg_value);


figure
ima_bin = imagesc(bin_lon_mat(1,:), bin_lat_mat(:,1), no2_bin);
set(ima_bin, 'AlphaData', ~isnan(no2_bin))
xlabel('distance (km)')
ylabel('distance (km)')
title('TROPOMI data gridded (filled)')
set(gca, 'YDir', 'normal','FontSize', 12,'fontweight','bold','FontName', 'Times New Roman')
pbaspect([4 1 1])
axis tight;

disp('Does the plume along the wind direction look east-west?')
cor_plu_dir_usr = input('Choose yes/no: ', 's');

if strcmp(cor_plu_dir_usr, 'no')
    error('You chose "no". Program terminated.');
end

plu_value = prctile(no2_bin(:), plu_per);

no2_bin_enh = no2_bin - plu_value;
no2_bin_enh(no2_bin_enh <= 0) = 0;
[row_plu, col_plu] = find (no2_bin_enh > 0);


l_con = bwlabeln(no2_bin_enh);
figure
imagesc(bin_lon_mat(1,:), bin_lat_mat(:,1), l_con)
set(ima_bin, 'AlphaData', ~isnan(l_con))
xlabel('distance (km)')
ylabel('distance (km)')
title('Plume classification')
set(gca, 'YDir', 'normal','FontSize', 12,'fontweight','bold','FontName', 'Times New Roman')
pbaspect([4 1 1])
colormap(jet(length(unique(l_con(~isnan(l_con))))));
colorbar;
axis tight;

l_con(l_con==0) = nan;
main_plu_l_con = mode(l_con, 'all');
l_con(isnan(l_con)) = 0;

disp(['The program find the plume ', num2str(main_plu_l_con), ' is the required source plume'])

cor_plu_usr = input('choose yes/no: ', 's');

if strcmp(cor_plu_usr, 'no')
    main_plu_l_con = input('choose the main plume based on Plume classification figure (1 or 2 or 3 or etc..): ');
end

for plu_xx = 1:length(row_plu)
    if l_con (row_plu(plu_xx), col_plu(plu_xx)) ~= main_plu_l_con
        no2_bin (row_plu(plu_xx), col_plu(plu_xx)) = bkg_value;
    elseif isnan(l_con (row_plu(plu_xx), col_plu(plu_xx)))
        no2_bin (row_plu(plu_xx), col_plu(plu_xx)) = bkg_value;
    end
end


figure
ima_bin_enh = imagesc(bin_lon_mat(1,:), bin_lat_mat(:,1), no2_bin);
set(ima_bin_enh, 'AlphaData', ~isnan(no2_bin))
xlabel('distance (km)')
ylabel('distance (km)')
title('TROPOMI data gridded (only source plume)')
set(gca, 'YDir', 'normal','FontSize', 12,'fontweight','bold','FontName', 'Times New Roman')
pbaspect([4 1 1])
axis tight;

end

