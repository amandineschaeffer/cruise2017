clear all 
close all

FILE = '../DATA/in2017_v04uwy5min.csv';
C=importdata(FILE,',',1);
% Get HEADER TEXT
header = C.textdata(1,:)';
% Get DAY/MONTH/YEAR
datecell = flipud(C.textdata(2:end,1)); 
DMY = cell2mat(datecell);
% Get HH/MM/SS
timecell = flipud(C.textdata(2:end,2));
HMS = cell2mat(timecell);
% Make TIME % dt = 5 min
for t=1:length(HMS)
   time0 = [char(DMY(t,:)) '-' char(HMS(t,:))] ;
   time(t)=datenum(time0,'dd-mmm-yyyy-HH:MM:SS');
end
datestr(time)
% Get VARIABLES [minus 2 since two less data than textdata]
heading = flipud(C.data(:,9));
lat = flipud(C.data(:,1)); % Lat
lon = flipud(C.data(:,2)); % Lon
% Salinity
sal = flipud(C.data(:,63-2)); sal(sal < 35.3) = NaN;
% Temperature water
temp = flipud(C.data(:,65-2));
% Port True Wind Speed (kts)
wsp_kts = flipud(C.data(:,27-2));
% Port True Wind Dir 
wdp = flipud(C.data(:,35-2));
% Starboard True Wind Speed (kts)
wss_kts = flipud(C.data(:,31-2));
% Starboard True Wind Dir 
wds = flipud(C.data(:,39-2));
% Get u/v from r/theta WRONG BECAUSE WRONG CONVENTION!
% [us,vs] = pol2cart(deg2rad(wds),wss);
% [up,vp] = pol2cart(deg2rad(wdp),wsp);

%% Convert to m/s
wsp = convvel(wsp_kts,'kts','m/s');
wss = convvel(wss_kts,'kts','m/s');

%% Note wind direction is in meteorological convention, i.e. wind blowing FROM the north is 0edg, then clockwise.
[us,vs]=compass2uv(wds,wss);
[up,vp]=compass2uv(wdp,wsp);
% Wind stress
[WIND_us_stress,WIND_vs_stress]=ra_windstr(us,vs);
[WIND_up_stress,WIND_vp_stress]=ra_windstr(up,vp);
% % Boat heading
% [u,v]=compass2uv(heading,1);
% [theta,rho] = cart2pol(u,v);
%%% Rotate velocities perpendicular from transect track
for i=2:length(time)
 [dist(i),phaseangle(i)] = sw_dist([lat(i),lat(i-1)],[lon(i),lon(i-1)],'km') ;
end
 rot_deg_angle=-(-phaseangle+90); % Note for CH it is -22deg base on NS axis;
if abs(rot_deg_angle)>90
    rot_deg_angle=rot_deg_angle+180;
end
WIND_up_stress_transect=cosd(rot_deg_angle).*WIND_up_stress'+sind(rot_deg_angle).*WIND_vp_stress'; % across-shelf
WIND_vp_stress_transect=-sind(rot_deg_angle).*WIND_up_stress'+cosd(rot_deg_angle).*WIND_vp_stress'; % along-shelf


%% Density
pres = 6+0*sal; %?????????????? 6 m depth???
% Absolute salinity
[SA, in_ocean] = gsw_SA_from_SP(sal,pres,lon,lat);
% Conservative temperature
CT = gsw_CT_from_t(sal,temp,pres);
% Density
dens=sw_dens(sal,temp,pres);    % Not with SA and CT... careful!
DENS=sw_dens(SA,CT,pres);

%% Distance
dist_grid_dx = deg2km(distance(lat(2:end),lon(2:end),lat(1:end-1),lon(1:end-1)));   % Careful, km!
dist_grid = [0 cumsum(dist_grid_dx)']';  % Careful small....
% buoyancy:
g = gsw_grav(lat,pres);
rho_0 = nanmean(dens);
buoyancy = - g.* (dens - rho_0)./ rho_0;
M2 = diff(buoyancy,1,1)./(dist_grid_dx*1000);  % in metres!
M2(dist_grid_dx<0.2)=NaN;   % in km!
M2(end+1) = M2(end); % return to original grid

% Ekman Buoyancy flux (D'asaro 2011)
f = gsw_f(lat);
EBF = WIND_up_stress_transect'.* M2./f./rho_0;  % or WIND_vp_stress: cross-transect is almost northward


%%
% Index Wind during Triaxus Tow
index_Triaxus = find(time > datenum('12-Sep-2017 14:38:50') & time < datenum('12-Sep-2017 22:33:10'));
% Index Wind 2 days before and after tow
index_Around = find(time > datenum('09-Sep-2017 14:38:50') & time < datenum('14-Sep-2017 22:33:10'));

% Plot Wind Components
figure;
subplot(2,1,1)
plot(time,us,time,up);datetick;axis tight
legend('Starboard U','Port U')
subplot(2,1,2)
plot(time,vs,time,vp);datetick;axis tight
legend('Starboard V','Port V')

% Plot T S 
figure; set(gcf, 'position',[2349   287    987    600],'PaperPositionMode','auto' , 'Color', 'w')
subplot(3,1,1)
plot(time(index_Around),temp(index_Around));hold on
plot(time(index_Triaxus),temp(index_Triaxus),'linewidth',3);datetick
title('Underway temperature'); ylabel('[^oC]')
subplot(3,1,2)
plot(time(index_Around),sal(index_Around));hold on
plot(time(index_Triaxus),sal(index_Triaxus),'linewidth',3);datetick
title('Underway salinity'); ylabel('[PSU]'); 
subplot(3,1,3)
plot(time(index_Around),EBF(index_Around));hold on
hold on
plot(time(index_Triaxus),EBF(index_Triaxus),'linewidth',3);datetick
title('Ekman buoyancy flux'); xlabel('Date')
export_fig(['../PLOTS/plot_underway_T_S_EBF.png'],'-opengl','-m2','-q101')

figure; set(gcf, 'position',[2349   287    987    300],'PaperPositionMode','auto' , 'Color', 'w')
plot(time(index_Around),EBF(index_Around));hold on
hold on
plot(time(index_Triaxus),EBF(index_Triaxus),'linewidth',3);datetick
ylabel(''); xlabel('Date')

figure; set(gcf, 'position',[2349   287    987    600],'PaperPositionMode','auto' , 'Color', 'w')
subplot(2,1,1)
plot(time(index_Around),us(index_Around),time(index_Around),up(index_Around));hold on
plot(time(index_Triaxus),us(index_Triaxus),time(index_Triaxus),up(index_Triaxus),'linewidth',3);datetick
legend('Starboard U','Port U'); ylabel('[m s^-^1]')
hold on;plot(time,zeros(1,length(time)),'k--')
title('Wind')
subplot(2,1,2)
plot(time(index_Around),vs(index_Around),time(index_Around),vp(index_Around));hold on
plot(time(index_Triaxus),vs(index_Triaxus),time(index_Triaxus),vp(index_Triaxus),'linewidth',3);datetick
legend('Starboard V','Port V'); ylabel('[m s^-^1]')
hold on;plot(time,zeros(1,length(time)),'k--')
export_fig(['../PLOTS/plot_underway_wind.png'],'-opengl','-m2','-q101')


% check wind direction: during the triaxus, almost only v component = across-transect
% figure;
% plot(time(index_Around),WIND_up_stress_transect(index_Around),time(index_Around),WIND_vp_stress(index_Around));hold on
% plot(time(index_Triaxus),WIND_up_stress_transect(index_Triaxus),time(index_Triaxus),WIND_vp_stress(index_Triaxus),'linewidth',3);datetick
% legend('Starboard U','Port U')
% hold on;plot(time,zeros(1,length(time)),'k--')

figure
quiver(lon(index_Around),lat(index_Around),us(index_Around),vs(index_Around),1,'g');hold on
quiver(lon(index_Around),lat(index_Around),up(index_Around),vp(index_Around),1,'r');hold on
quiver(lon(index_Triaxus),lat(index_Triaxus),us(index_Triaxus),vs(index_Triaxus),1,'g','LineWidth',2);hold on
quiver(lon(index_Triaxus),lat(index_Triaxus),up(index_Triaxus),vp(index_Triaxus),1,'r','LineWidth',2);hold on
Y=get(gca,'ylim');
set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
box on

%% save
WIND = struct('time',time,'us',us,'vs',vs,'up',up,'vp',vp,'lat',lat,'lon',lon);
save save_wind_raw.mat WIND