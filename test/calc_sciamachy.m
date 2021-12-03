%% 
%%% ======================================================================
%%% calc_sciamachy.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% collect observations for the current model timestep
s = 1;
for k=1:3:length(t)
    if t_step(nstep)<t(k) && t(k)<t_step(nstep+1)
        lat1(s)   = lat(k);
        lon1(s)   = lon(k);
        hm1(s)    =  hm(k);
        ch41(s)   = ch4(k);
        sig1(s)   = sig(k);
        avk1(:,s) = avk(:,k);
        plv1(:,s) = plv(:,k);
        pwe1(:,s) = pwe(:,k);
        fgs1(:,s) = fgs(:,k);
        s=s+1;
    end
end

if exist('lat1','var') % if observations are not available, no need to process

%%% index to limit data to Northern hemisphere
    in = ingeoquad(lat1,lon1,latlim,lonlim);
    lat_sat_N = lat1(in);
    lon_sat_N = lon1(in);
    sig_sat_N = sig1(:,in);
    ch4_sat_N = ch41(:,in);
    plv_sat_N = plv1(:,in);
    pwe_sat_N = pwe1(:,in);
    fgs_sat_N = fgs1(:,in);
    avk_sat_N = avk1(:,in);

    SCIAMACHY_TIME = [t_step(nstep), t_step(nstep+1)]
    if ~isempty(lat_sat_N)
    	disp('QC SCIAMACHY processing ...')
        
%%% quality control to remove the outlier whose departure from mean > 3std. 
        [indx_out,L_out,U_out,C_out]    = isoutlier(ch4_sat_N,'mean','ThresholdFactor',3);
        ch4_sat1  = ch4_sat_N(~indx_out);
        sig_sat1  = sig_sat_N(~indx_out);
        lat_sat1  = lat_sat_N(~indx_out);
        lon_sat1  = lon_sat_N(~indx_out);
        plv_sat1  = plv_sat_N(:,~indx_out);
        pwe_sat1  = pwe_sat_N(:,~indx_out);
        fgs_sat1  = fgs_sat_N(:,~indx_out);
        avk_sat1  = avk_sat_N(:,~indx_out);
        [xproj_sat1,yproj_sat1] = ll2psn(lat_sat1,lon_sat1,'TrueLat',45,'EarthRadius',6370000,...
    'Eccentricity',1e-30,'meridian',-98);	
	kx_indx  = find(min(xproj_hcmaq(:)) < xproj_sat1 & xproj_sat1 < max(xproj_hcmaq(:)));
        ky_indx  = find(min(yproj_hcmaq(:)) < yproj_sat1 & yproj_sat1 < max(yproj_hcmaq(:)));
        kxy_indx = intersect(kx_indx,ky_indx);
        ch4_sat  = ch4_sat1(kxy_indx);
        f_rep    = 0;
        sig_sat  = sig_sat1(kxy_indx) + f_rep * (0.015 * ch4_sat);
        %sig_sat  = sig_sat1(kxy_indx);
        lat_sat  = lat_sat1(kxy_indx);
        lon_sat  = lon_sat1(kxy_indx);
        plv_sat  = plv_sat1(:,kxy_indx);
        pwe_sat  = pwe_sat1(:,kxy_indx);
        fgs_sat  = fgs_sat1(:,kxy_indx);
        avk_sat  = avk_sat1(:,kxy_indx);
        
%%% bias correction (latest update)
    	lat_sat_rad = deg2rad(lat_sat);
    	lon_sat_rad = deg2rad(lon_sat);
    	[x_sat,y_sat,z_sat] = sph2cart(lon_sat_rad,lat_sat_rad,r);
    	r_sat  = [x_sat(:),y_sat(:),z_sat(:)];
    	[xproj_sat,yproj_sat] = ll2psn(lat_sat,lon_sat,'TrueLat',45,'EarthRadius',6370000,...
    	'Eccentricity',1e-30,'meridian',-98);
    
%%% if no observations exist make variables empty 
    else
        disp('SCIAMACHY Unavailable')
        lat_sat = [];
        lon_sat = [];
        ch4_sat = [];
        sig_sat = [];
        plv_sat = [];
        pwe_sat = [];
        fgs_sat = [];
        avk_sat = [];
    end
%%% if no observations exist make variables empty     
else
disp('SCIAMACHY Unavailable')
lat_sat = [];
lon_sat = [];
ch4_sat = [];
sig_sat = [];
plv_sat = [];
pwe_sat = [];
fgs_sat = [];
avk_sat = [];
end
%%% =======================================================================
%%% END
%%% =======================================================================