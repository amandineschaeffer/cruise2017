% Plot Triaxus
clear;clc
close all

path(path,genpath('/home/amandine/Logiciels/matlab/seawater'));
path(path,genpath('/home/amandine/Logiciels/matlab/GSW_V305_toolbox'));
path(path,genpath('/home/amandine/Logiciels/matlab/altmany-export_fig-7720793'));

% TRIAXUS__________________________________________________________________
% addpath ./matlab_scripts/Triaxus/
[lon,lat,time,pres,temp,sal,oxy,chl,tur,cdom] = process_TRIAXUS;
depth = - sw_dpth(pres,lat);

figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(4,1,1)
scatter(lon,depth,40,temp,'filled');
set(gca,'fontsize',12)
colormap jet
ht=colorbar; title(ht,'Temperature (^o C)')
ylabel('Depth (m)')

subplot(4,1,2)
scatter(lon,depth,40,sal,'filled');
set(gca,'fontsize',12)
colormap jet
hs=colorbar; title(hs,'Salinity (PSU)')
ylabel('Depth (m)')

subplot(4,1,3)
scatter(lon,depth,40,oxy,'filled');
set(gca,'fontsize',12)
colormap jet
ho=colorbar; title(ho,'Oxygen (uM/l)')
ylabel('Depth (m)')

subplot(4,1,4)
scatter(lon,depth,40,chl,'filled');
set(gca,'fontsize',12)
colormap jet
hc=colorbar; title(hc,'Chlorophyll (kg m^{-3}?)')
ylabel('Depth (m)')
xlabel('Longitude')
export_fig(['../PLOTS/plot_triaxus_RAW.png'],'-opengl','-m2','-q101')

% Grid over pres or depth (change sign)
dist = deg2km(distance(lon,lat,lon(1),lat(1)));
ind_bottom = find(pres(3:end)<pres(2:end-1) & pres(1:end-2)< pres(2:end-1) & pres(2:end-1) > 50)-1; % bottom >140m
ind_surf = find(pres(3:end)>pres(2:end-1) & pres(1:end-2)>pres(2:end-1) )+1 ; %surface

figure;
plot(lon,depth)
hold on
scatter(lon(ind_bottom),depth(ind_bottom))
scatter(lon(ind_surf),depth(ind_surf))
% scatter(time(ind_bottom),depth(ind_bottom))
% scatter(time(ind_surf),depth(ind_surf))

depth_max = 340
% Average each down / up cast
if length(ind_bottom) ~= length(ind_bottom), disp('Warning number full up / down cast different'); end
for i=1:length(ind_bottom)-1
    ind_cast = i
    depth_grid = [-depth_max:1];   % New grid
    % Average time, lat ,lon
    lat_grid(i) = mean(lat(ind_bottom(i):ind_bottom(i+1)));
    lon_grid(i) = mean(lon(ind_bottom(i):ind_bottom(i+1)));
    time_grid(i) = mean(time(ind_bottom(i):ind_bottom(i+1)));
    
    % DOWN cast interp;
    ind_down = ind_bottom(i):ind_surf(i);
    depth_ori = round(depth(ind_down));
    temp_ori = temp(ind_down);
    sal_ori = sal(ind_down);
    oxy_ori = oxy(ind_down);
    chl_ori = chl(ind_down);
    tur_ori = tur(ind_down);
    cdom_ori = cdom(ind_down);
    temp_grid1 = fun_interpolate_cast(-depth_ori,temp_ori,-depth_grid); % function amandine
    sal_grid1 = fun_interpolate_cast(-depth_ori,sal_ori,-depth_grid); % function amandine
    oxy_grid1 = fun_interpolate_cast(-depth_ori,oxy_ori,-depth_grid); % function amandine
    chl_grid1 = fun_interpolate_cast(-depth_ori,chl_ori,-depth_grid); % function amandine
    tur_grid1 = fun_interpolate_cast(-depth_ori,tur_ori,-depth_grid); % function amandine
    cdom_grid1 = fun_interpolate_cast(-depth_ori,cdom_ori,-depth_grid); % function amandine
    
    % UP cast interp;
    ind_up = ind_surf(i):ind_bottom(i+1);
    depth_ori = round(depth(ind_up));
    temp_ori = temp(ind_up);
    sal_ori = sal(ind_up);
    oxy_ori = oxy(ind_up);
    chl_ori = chl(ind_up);
    tur_ori = tur(ind_up);
    cdom_ori = cdom(ind_up);
    temp_grid2 = fun_interpolate_cast(-depth_ori,temp_ori,-depth_grid); % function amandine
    sal_grid2 = fun_interpolate_cast(-depth_ori,sal_ori,-depth_grid); % function amandine
    oxy_grid2 = fun_interpolate_cast(-depth_ori,oxy_ori,-depth_grid); % function amandine
    chl_grid2 = fun_interpolate_cast(-depth_ori,chl_ori,-depth_grid); % function amandine
    tur_grid2 = fun_interpolate_cast(-depth_ori,tur_ori,-depth_grid); % function amandine
    cdom_grid2 = fun_interpolate_cast(-depth_ori,cdom_ori,-depth_grid); % function amandine
    
    % Average:
    TEMP_grid(i,:) = nanmean([temp_grid1',temp_grid2'],2);
    SAL_grid(i,:) = nanmean([sal_grid1',sal_grid2'],2);
    OXY_grid(i,:) = nanmean([oxy_grid1',oxy_grid2'],2);
    CHL_grid(i,:) = nanmean([chl_grid1',chl_grid2'],2);
    TUR_grid(i,:) = nanmean([tur_grid1',tur_grid2'],2);
    CDOM_grid(i,:) = nanmean([cdom_grid1',cdom_grid2'],2);
    %   figure; plot(temp(ind_bottom(i):ind_bottom(i+1)),depth(ind_bottom(i):ind_bottom(i+1)));    hold on ;   plot(TEMP_grid(i,:),depth_grid) % Check
    
end
% MORE
dist_grid_dx = deg2km(distance(lat_grid(2:end),lon_grid(2:end),lat_grid(1:end-1),lon_grid(1:end-1)));   % Careful, km!
dist_grid = [0 cumsum(dist_grid_dx)];   % Careful, km!
[DEPTH_grid,CAST_grid] = meshgrid(depth_grid,[1:length(ind_surf)-1]);
[~,LAT_grid] = meshgrid(depth_grid,lat_grid);
[~,LON_grid] = meshgrid(depth_grid,lon_grid);
[~,TIME_grid] = meshgrid(depth_grid,time_grid);
[~,DIST_grid] = meshgrid(depth_grid,dist_grid);
PRES_grid = sw_pres(-DEPTH_grid,ones(size(LAT_grid)).*nanmean(nanmean(LAT_grid,1),2));

% Absolute salinity
[SA_grid, in_ocean] = gsw_SA_from_SP(SAL_grid,PRES_grid,LON_grid,LAT_grid);
% Conservative temperature
CT_grid = gsw_CT_from_t(SAL_grid,TEMP_grid,PRES_grid);
% Density
dens=sw_dens(sal,temp,pres);    % Not with SA and CT... careful!
DENS_grid=sw_dens(SA_grid,CT_grid,PRES_grid);


% TEMP
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,temp,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Temperature (^oC)')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,TEMP_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% OXY
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,oxy,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Oxygen (uM/l)')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,OXY_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% SAL
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,sal,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Salinity')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,SAL_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% TURB
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,tur,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Turbidity')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,TUR_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% CHL
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,chl,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Chlorophyll (kg m^{-3}?)')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,CHL_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% CDOM
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,cdom,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('CDOM')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,CDOM_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])
% DENS
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
scatter(time,depth,40,dens,'filled'); colormap(jet(10)) ; ht=colorbar; ylabel(ht,'raw'); ylabel('Depth (m)')
datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
title('Density')
subplot(2,1,2,'align')
contourf(TIME_grid,DEPTH_grid,DENS_grid,10,'k','linewidth',2);
colormap(jet(10)) ; ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)');xlabel('Distance (km)')
datetick('x');  xlim([time(end) time(1)]);  ylim([-350 0])



%% N2: buoyancy (Brunt-Vaisala) frequency squared (N^2)
% One way...
[N2_aaa, p_mid_aaa] = gsw_Nsquared(SA_grid',CT_grid',PRES_grid',LAT_grid'); % at mid-points pressure
N2_grid = interp2(p_mid_aaa',DIST_grid(:,1:end-1),N2_aaa',PRES_grid,DIST_grid); % interpolate
% other way...
[N2_v2_aaa,q,p_ave_aaa] = sw_bfrq(SAL_grid',TEMP_grid',PRES_grid',LAT_grid');  % at mid-points pressure
N2_v2_grid = interp2(p_ave_aaa',DIST_grid(:,1:end-1),N2_v2_aaa',PRES_grid,DIST_grid); % interpolate
% PV_v2_grid = interp2(p_ave_aaa',DIST_grid(:,1:end-1),q',PRES_grid,DIST_grid); % interpolate

% N2
figure; set(gcf, 'position',[2349   287    987    600],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1);
pcolor(DIST_grid,DEPTH_grid,N2_grid);    shading flat;  title(['Brunt-Vaisala frequency squared (N^2)']);
colormap(jet); colorbar; ylabel('Depth (m)'); ; xlabel('Distance (km)')
% caxis([-1e-4 1e-4])
subplot(2,1,2);
pcolor(DIST_grid,DEPTH_grid,N2_grid);    shading flat;  title(['Brunt-Vaisala frequency squared (N^2)']);
colormap(jet); colorbar;  ylabel('Depth (m)');; xlabel('Distance (km)')
export_fig(['../PLOTS/plot_N2.png'],'-opengl','-m2','-q101')


%% M2:
[M2_aaa,p_ave_aaa] = AM_sw_M2(SAL_grid',TEMP_grid',PRES_grid',LAT_grid',LON_grid');
M2_grid = interp2(p_ave_aaa',DIST_grid(1:end-1,:),M2_aaa',PRES_grid,DIST_grid); % interpolate
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,M2_grid);    shading flat;  title(['Horizontal buoyancy gradient (M^2)']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
% caxis([-1e-6 1e-6])
export_fig(['../PLOTS/plot_M2.png'],'-opengl','-m2','-q101')


%% Balanced Richardson nb:f2 N2 /(M4)
% Coriolis

f = gsw_f(LAT_grid);
Ri_b = (f.^2 .* N2_grid) ./ M2_grid.^2;
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,Ri_b);    shading flat;  title([' Balanced Richardson nb']);
colormap(jet); colorbar; ylabel('Depth (m)'); xlabel('Distance (km)')
caxis([0 1])
export_fig(['../PLOTS/plot_Ri_balanced.png'],'-opengl','-m2','-q101')

Phi_RI_b = atan(-1./Ri_b);   % radians
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,Phi_RI_b*180/pi);    shading flat;  title(['\phi Ri_b (degrees)']);
colormap(jet); colorbar; datetick('x');ylabel('Depth (m)'); xlabel('Distance (km)')
export_fig(['../PLOTS/plot_Phi_Rb.png'],'-opengl','-m2','-q101')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET ADCP

% % SHIPBOARD ADCP___________________________________________________________
% addpath ./matlab_scripts/ADCP/
[lon_adcp,lat_adcp,depth_adcp,time_adcp,u_adcp,v_adcp,lon_track_adcp,lat_track_adcp,u_track_adcp,v_track_adcp,div_adcp,vort_adcp] = process_ADCP();
% Add time u, v


% Interpolate into Triaxus grid
figure; scatter(LON_grid(:,1),LAT_grid(:,1)); hold on; scatter(lon_adcp,lat_adcp)
dist_grid_dx_adcp = deg2km(distance(lat_adcp(2:end),lon_adcp(2:end),lat_adcp(1:end-1),lon_adcp(1:end-1)));
dist_grid_adcp = [0 cumsum(dist_grid_dx_adcp)];
dist_adcp_rel_triaxus = deg2km(distance(lat_adcp(1),lon_adcp(1),LAT_grid(1,1),LON_grid(1,1)));
dist_grid_adcp_rel2Triaxus = dist_grid_adcp + dist_adcp_rel_triaxus;
% Track
dist_grid_dx_track_adcp = deg2km(distance(lat_track_adcp(2:end),lon_track_adcp(2:end),lat_track_adcp(1:end-1),lon_track_adcp(1:end-1)));
dist_grid_track_adcp = [0 cumsum(dist_grid_dx_track_adcp)];
dist_track_adcp_rel_triaxus = deg2km(distance(lat_track_adcp(1),lon_track_adcp(1),LAT_grid(1,1),LON_grid(1,1)));
dist_grid_track_adcp_rel2Triaxus = dist_grid_track_adcp + dist_track_adcp_rel_triaxus;


% dist_ADCP = deg2km(distance(lon_adcp,lat_adcp,lon_grid(1),lat_grid(1)));
% dist_vorti_ADCP = deg2km(distance(lon_vorti_adcp,lat_vorti_adcp,lon_grid(1),lat_grid(1)));
% dist_track_ADCP = deg2km(distance(lon_track_adcp,lat_track_adcp,lon_grid(1),lat_grid(1)));

% U_ADCP_grid = interp2(dist_ADCP,pres_adcp,u_adcp,DIST_grid,PRES_grid);
% V_ADCP_grid = interp2(dist_ADCP,pres_adcp,v_adcp,DIST_grid,PRES_grid);
U_track_ADCP_grid = interp2(dist_grid_adcp_rel2Triaxus,depth_adcp,u_track_adcp,DIST_grid,DEPTH_grid);
V_track_ADCP_grid = interp2(dist_grid_adcp_rel2Triaxus,depth_adcp,v_track_adcp,DIST_grid,DEPTH_grid);
VORTI_ADCP_grid = interp2(dist_grid_track_adcp_rel2Triaxus,depth_adcp,vort_adcp,DIST_grid,DEPTH_grid);
DIV_ADCP_grid = interp2(dist_grid_track_adcp_rel2Triaxus,depth_adcp,div_adcp,DIST_grid,DEPTH_grid);

% CHECK
% % U_ADCP
% figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
% subplot(2,1,1)
% pcolor(lon_adcp,depth_adcp,u_adcp) ; shading flat;  ht=colorbar; ylabel(ht,'ori'); ylabel('Depth (m)')
% datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])
% title('U adcp')
% subplot(2,1,2,'align')
% pcolor(TIME_grid,DEPTH_grid,U_ADCP_grid) ; shading flat;  ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)')
% datetick('x');  xlim([time(end) time(1)]);   ylim([-350 0])

% U_ADCP ACROSS
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
pcolor(lon_adcp,depth_adcp,v_track_adcp) ; shading flat;  ht=colorbar; ylabel(ht,'ori'); ylabel('Depth (m)')
ylim([-350 0])
title('U across adcp')
subplot(2,1,2,'align')
pcolor(LON_grid,DEPTH_grid,V_track_ADCP_grid) ; shading flat;  ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)')
ylim([-350 0])

% VORTI ACROSS
figure; set(gcf, 'position',[2349   287    987    646],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
pcolor(lon_track_adcp,depth_adcp,vort_adcp) ; shading flat;  ht=colorbar; ylabel(ht,'ori'); ylabel('Depth (m)')
ylim([-350 0])
title('Vorti adcp')
subplot(2,1,2,'align')
pcolor(LON_grid,DEPTH_grid,VORTI_ADCP_grid) ; shading flat;  ht=colorbar; ylabel(ht,'interpolated'); ylabel('Depth (m)')
ylim([-350 0])

% VORTI ACROSS
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,VORTI_ADCP_grid./abs(f)) ; shading flat;  ht=colorbar; ylabel(ht,'vorticity / |f|'); ylabel('Depth (m)')
caxis([-1.5 1.5]); colormap(r_b); xlabel('Distance (km)')
hold on
contour(DIST_grid,DEPTH_grid,TEMP_grid,30,'k','linewidth',2)
ylim([-350 0])
export_fig(['../PLOTS/plot_vorti_tempCONTOURS.png'],'-opengl','-m2','-q101')



%% Richardson nb: N2 / (du/dz)^2
% [aaa,DUDz] = gradient(U_track_ADCP_grid,1,-1) ;  %  DUDz = diff(U_across_ADCP_grid,1,2) same %[FX,FY] = GRADIENT(F) returns the numerical gradient of the
% %   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) direction.
% [aaa,DVDz] = gradient(V_track_ADCP_grid,1,-1) ;
DUDz = diff(U_track_ADCP_grid,1,2)./diff(DEPTH_grid,1,2); 
DVDz = diff(V_track_ADCP_grid,1,2)./diff(DEPTH_grid,1,2);
grid_aaa = (DEPTH_grid(:,2:end)+DEPTH_grid(:,1:end-1))/2;
DUDz_grid = interp2(DIST_grid(:,1:end-1)',grid_aaa',DUDz',DIST_grid,DEPTH_grid);
DVDz_grid = interp2(DIST_grid(:,1:end-1)',grid_aaa',DVDz',DIST_grid,DEPTH_grid);
Ri=N2_grid./(DUDz_grid.^2 + DVDz_grid.^2);
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,Ri); shading flat;  ht=colorbar; ylabel(ht,''); ylabel('Depth (m)')
caxis([0 1]); xlabel('Distance (km)')
title(['Richardson Number (N^2 / (du/dz)^2)']); ylim([-350 0])
export_fig(['../PLOTS/plot_Ri.png'],'-opengl','-m2','-q101')




%%% Critical angle:
vorti_a = f + VORTI_ADCP_grid;
Ri_C = atan( - vorti_a./ f);
% %             Ri_C (f.* vorti_a <0) = NaN;

% Potential vorticity:
PV = (f + VORTI_ADCP_grid).* N2_grid - M2_grid.^2./ f;
Qvert = (f + VORTI_ADCP_grid).* N2_grid;
Qbc = - M2_grid.^2./ f;
Vorti_abs = (f + VORTI_ADCP_grid);

%
%              figure;
%              if i ==2
%                 set(gcf, 'position',[ 825    68   842   266],'PaperPositionMode','auto' , 'Color', 'w')
%              else
%                 set(gcf, 'position',[ 1229          68         438         266],'PaperPositionMode','auto' , 'Color', 'w')
%              end
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,PV*1E9); shading interp; title(['Potential Vorticity x 10^9 (s^-^3)']);
colormap(r_b); colorbar;  ylabel('Depth (m)');caxis([-5 5 ]);  xlabel('Distance (km)')
hold on ; [C,h] = contour(PV*1E9,[0 0],'k'); ylim([-350 0])
export_fig(['../PLOTS/plot_PV.png'],'-opengl','-m2','-q101')
%
figure; set(gcf, 'position',[2349   287    987    900],'PaperPositionMode','auto' , 'Color', 'w')
subplot(3,1,1)
pcolor(DIST_grid,DEPTH_grid,PV*1E9); shading interp; title(['Potential Vorticity x 10^9 (s^-^3)']);
colormap(r_b); colorbar;  ylabel('Depth (m)');caxis([-5 5 ]);  xlabel('Distance (km)')
hold on ; [C,h] = contour(PV*1E9,[0 0],'k'); ylim([-350 0])
subplot(3,1,2)
pcolor(DIST_grid,DEPTH_grid,Qvert); shading interp; title(['Qvert']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
% export_fig(['../PLOTS/plot_PV.png'],'-opengl','-m2','-q101')
subplot(3,1,3)
pcolor(DIST_grid,DEPTH_grid,Qbc); shading interp; title(['Qbc']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
 export_fig(['../PLOTS/plot_PV_terms.png'],'-opengl','-m2','-q101')
% export_fig(['../PLOTS/plot_PV.png'],'-opengl','-m2','-q101')
figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,Vorti_abs); shading interp; title(['Vorti_abs']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
% % export_fig(['../PLOTS/plot_PV.png'],'-opengl','-m2','-q101')
% figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
% pcolor(DIST_grid,DEPTH_grid,VORTI_ADCP_grid.*abs(f')); shading interp; title(['Vorti_abs']);
% colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)'); caxis([-10*1e-5 10*1e-5])
% % export_fig(['../PLOTS/plot_PV.png'],'-opengl','-m2','-q101')



figure; set(gcf, 'position',[2349   287    987    600],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1); pcolor(DIST_grid,DEPTH_grid,Phi_RI_b); shading flat; title(['\phi Ri_b']); caxis([-1.5 1.5])
colormap(jet(4)); colorbar; ylabel('Depth (m)'); 
subplot(2,1,2); plot(DIST_grid(:,1),Ri_C,'x')

INSTAB = Phi_RI_b;
INSTAB_SI = NaN.* Phi_RI_b;
INSTAB_cat = NaN.* Phi_RI_b;
ANTICYL = (VORTI_ADCP_grid./f<1);
CYCL = (VORTI_ADCP_grid./f>1);
% STABLE:
INSTAB_cat(Phi_RI_b > Ri_C) = 200;    % Stable
% SI / intertial
INSTAB_cat(ANTICYL & Phi_RI_b > - pi/4 & Phi_RI_b < Ri_C) = 100;    % mixed Inertial / SI
% SI:
INSTAB_cat(ANTICYL & Phi_RI_b > - pi/2 & Phi_RI_b < -pi/4) = 0;    % SI
INSTAB_cat(CYCL & Ri_C <-pi/4 & Phi_RI_b > - pi/2 & Phi_RI_b < Ri_C) = 0;    % SI
% SI / G:
INSTAB_cat(Phi_RI_b > - 3*pi/4 & Phi_RI_b < -pi/2) = -100;    % mixed GI / SI
% Gravitationnal:
INSTAB_cat(Phi_RI_b <- 3*pi/4) = -200;   % G

% for j=1:size(INSTAB_cat,1)
%     aaa = Phi_RI_b(j,:);
%     INSTAB_cat(j,aaa(2:end) > 0 & aaa(2:end) < Ri_C(j,:)) = 100;     %% ??
%     INSTAB_SI(j,aaa(2:end) > 0 & aaa(2:end) < Ri_C) = 1;     %% ??
%     %                 INSTAB_cat(j,aaa(1:end-1) > 0 & aaa(1:end-1) < Ri_C) = 100;
%     INSTAB_cat(j,aaa(2:end) > Ri_C) = NaN;     %% ??
%     INSTAB(j,aaa(2:end) > Ri_C) = NaN;
% end

figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,INSTAB_cat);colormap(jet(3)); shading flat; title(['Sub-mesoscale instabilities']);
caxis([-200 200])
colormap(jet(5)); colorbar('Ticks',[-200,-100,0,100,200], 'TickLabels',{'G','G/SI','SI','I/SI','Stable'})
ylabel('Depth (m)');
export_fig(['../PLOTS/plot_Instab_cat.png'],'-opengl','-m2','-q101')

%%% Frontogenesis
%-2b db/dx (du/dx * db/dx)
DUDx = diff(U_track_ADCP_grid,1,1)./diff(DIST_grid,1,1); 
grid_aaa = (DIST_grid(2:end,:)+DIST_grid(1:end-1,:))/2;
DUDx_grid = interp2(grid_aaa',DEPTH_grid(1:end-1,:)',DUDx',DIST_grid,DEPTH_grid);

HCL = -2 * M2_grid.* (DUDx_grid.* M2_grid) ;

figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,HCL); shading interp; title(['Frontogenesis: HCL term (horizontal strain and shearing)']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
caxis([-2e-14 2e-14])
export_fig(['../PLOTS/plot_HCL.png'],'-opengl','-m2','-q101')



%%% WIND
load save_wind_raw.mat
% Interp
WIND_up_interp1 = interp1(WIND.time,WIND.up,TIME_grid(:,1));
% figure; plot(WIND.time,WIND.up); hold on ; plot(TIME_grid(:,1), WIND_up_grid);
WIND_vp_interp1 = interp1(WIND.time,WIND.vp,TIME_grid(:,1));
WIND_us_interp1 = interp1(WIND.time,WIND.us,TIME_grid(:,1));
WIND_vs_interp1 = interp1(WIND.time,WIND.vs,TIME_grid(:,1));

[aaa,WIND_up_grid] = meshgrid(depth_grid,WIND_up_interp1);
[aaa,WIND_vp_grid] = meshgrid(depth_grid,WIND_vp_interp1);
[aaa,WIND_us_grid] = meshgrid(depth_grid,WIND_us_interp1);
[aaa,WIND_vs_grid] = meshgrid(depth_grid,WIND_vs_interp1);

[WIND_up_stress_grid,WIND_vp_stress_grid]=ra_windstr(WIND_up_grid,WIND_vp_grid);

% Ekman equivalent heat flux (Thompson 2016)
g = gsw_grav(LAT_grid,DEPTH_grid);
Cp = sw_cp(SAL_grid, TEMP_grid, PRES_grid); % Specific Heat Capacity  [J kg^-1 C^-1] 
alpha = sw_alpha(SAL_grid, TEMP_grid, PRES_grid); % Thermal expansion coeff (alpha) [degree_C.^-1]
Qek = - Cp.* M2_grid.* WIND_vp_stress_grid./ f./g ./ alpha;

figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,Qek); shading interp; title(['Qek']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
% caxis([-2e-14 2e-14])


% Ekman Buoyancy flux (D'asaro 2011)
rho_0 = nanmean(nanmean(DENS_grid));
EBF = WIND_vp_stress_grid.* M2_grid./f./rho_0;  % or WIND_vs_interp1 or get right angle...

figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
pcolor(DIST_grid,DEPTH_grid,EBF); shading interp; title(['EBF']);
colormap(r_b); colorbar;  ylabel('Depth (m)'); xlabel('Distance (km)')
% caxis([-2e-14 2e-14])




