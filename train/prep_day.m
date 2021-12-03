%% 
%%% ======================================================================
%%% prep_day.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% read lat and lon from H-CMAQ grid
hcmaq_lat  = double(ncread(Lname{i},'LAT')); % Lname --> read_filename_hemi
hcmaq_lon  = double(ncread(Lname{i},'LON')); % Lname --> read_filename_hemi
latlim  = [1 90];lonlim  = [-180 180]; % range of Northern hemisphere

%%% read surface and levels pressure and sigma levels HCMAQ 
hcmaq_p_s_25  = double(ncread(Mname_cro2d{i},'PRSFC')); % surface pressure 
hcmaq_siglvl  = ncreadatt(Mname_cro3d{i},'/','VGLVLS'); % sigma level value

%%% X-Y grided H-CMAQ on projection plane (ll2psn.m should be in current directory)
[x1,y1] = ll2psn(hcmaq_lat,hcmaq_lon,'TrueLat',45,'EarthRadius',6370000,...
    'Eccentricity',1e-30,'meridian',-98); % lat/lon to polarstereographic plane 
x1 = round(round(x1,6,'significant')); % rounding to ensure a regualr grid
y1 = round(round(y1,6,'significant')); % rounding to ensure a regualr grid
[xproj_hcmaq,yproj_hcmaq] = ndgrid(x1(1,1):+108000:x1(187,187),y1(1,1):...
    +108000:x1(187,187)); % create regualr grids of H-CMAQ from x1 and y1

%%% X-Y-Z from spherical to cartesian coordinate HCMAQ
hcmaq_lat_rad = deg2rad(hcmaq_lat);
hcmaq_lon_rad = deg2rad(hcmaq_lon);
r             = 6370000; % radius of Earth given in CMAQ as 6370000 m
[x_hcmaq,y_hcmaq,z_hcmaq] = sph2cart(hcmaq_lon_rad,hcmaq_lat_rad,r); %
r_hcmaq = [x_hcmaq(:),y_hcmaq(:),z_hcmaq(:)];

%%% =======================================================================
%%% END
%%% =======================================================================