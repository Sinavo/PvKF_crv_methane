%%
%%% ======================================================================
%%% calc_gosat_crossv.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

%%% read satellite observations from bias corrected file (.mat)
lat_sat  =  nonzeros(lat_sat_all_train(:,step_size*(i-1)+nstep));
lon_sat  =  nonzeros(lon_sat_all_train(:,step_size*(i-1)+nstep));
ch4_sat  =  ch4_sat_all_train(1:length(lat_sat),step_size*(i-1)+nstep);
sig_sat1 =  sig_sat_all_train(1:length(lat_sat),step_size*(i-1)+nstep);
plv_sat  =  plv_sat_all_train(:,1:length(lat_sat),step_size*(i-1)+nstep);
pwe_sat  =  pwe_sat_all_train(:,1:length(lat_sat),step_size*(i-1)+nstep);
fgs_sat  =  fgs_sat_all_train(:,1:length(lat_sat),step_size*(i-1)+nstep);
avk_sat  =  avk_sat_all_train(:,1:length(lat_sat),step_size*(i-1)+nstep);
%f_rep   = 3;
%sig_sat = sig_sat1 + f_rep * (0.007 * ch4_sat); % addetive obs. error
f_rep    = 0.45;
sig_sat  = f_rep.*(sig_sat1); % multiplicitive obs. error

%%% X-Y grided H-CMAQ on projection plane (ll2psn.m should be in current directory)
lat_sat_rad = deg2rad(lat_sat);
lon_sat_rad = deg2rad(lon_sat);
[x_sat,y_sat,z_sat] = sph2cart(lon_sat_rad,lat_sat_rad,r);
r_sat  = [x_sat(:),y_sat(:),z_sat(:)];
[xproj_sat,yproj_sat] = ll2psn(lat_sat,lon_sat,'TrueLat',45,'EarthRadius',6370000,...
    'Eccentricity',1e-30,'meridian',-98);
%%% =======================================================================
%%% END
%%% =======================================================================
