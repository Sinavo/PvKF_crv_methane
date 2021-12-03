%% 
%%% ======================================================================
%%% calc_gosat_offline.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% collect observations for the current model timestep
s = 1;
for k=1:1:length(t)
    if t_step(nstep)<t(k) && t(k)<t_step(nstep+1)
        lat1(s)   = lat(k);
        lon1(s)   = lon(k);
        tim1(s)   = tim(k);
        sen1(s)   = sen(k);
        sol1(s)   = sol(k);
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
    tim_sat_N = tim1(in);
    sen_sat_N = sen1(in);
    sol_sat_N = sol1(in);
    sig_sat_N = sig1(:,in);
    ch4_sat_N = ch41(:,in);
    plv_sat_N = plv1(:,in);
    pwe_sat_N = pwe1(:,in);
    fgs_sat_N = fgs1(:,in);
    avk_sat_N = avk1(:,in);

    GOSAT_TIME = [t_step(nstep), t_step(nstep+1)]
    if ~isempty(lat_sat_N)
    	disp('QC GOSAT processing ...')
        
%%% quality control to remove the outlier whose departure from mean > 3std. 
        [indx_out,L_out,U_out,C_out]    = isoutlier(ch4_sat_N,'mean','ThresholdFactor',3);
        ch4_sat1  = ch4_sat_N(~indx_out);
        sig_sat1  = sig_sat_N(~indx_out);
        lat_sat1  = lat_sat_N(~indx_out);
        lon_sat1  = lon_sat_N(~indx_out);
        tim_sat1  = tim_sat_N(~indx_out);
        sen_sat1  = sen_sat_N(~indx_out);
        sol_sat1  = sol_sat_N(~indx_out);
        plv_sat1  = plv_sat_N(:,~indx_out);
        pwe_sat1  = pwe_sat_N(:,~indx_out);
        fgs_sat1  = fgs_sat_N(:,~indx_out);
        avk_sat1  = avk_sat_N(:,~indx_out);
        [xproj_sat1,yproj_sat1] = ll2psn(lat_sat1,lon_sat1,'TrueLat',45,'EarthRadius',6370000,...
    'Eccentricity',1e-30,'meridian',-98);	
	kx_indx  = find(min(xproj_hcmaq(:)) < xproj_sat1 & xproj_sat1 < max(xproj_hcmaq(:)));
        ky_indx  = find(min(yproj_hcmaq(:)) < yproj_sat1 & yproj_sat1 < max(yproj_hcmaq(:)));
        kxy_indx = intersect(kx_indx,ky_indx);
        ch4_sat2  = ch4_sat1(kxy_indx);
        %sig_sat = sig_sat1(kxy_indx);
        f_rep    = 0;
        sig_sat2  = sig_sat1(kxy_indx) + f_rep * (0.007 * ch4_sat2);
        lat_sat2  = lat_sat1(kxy_indx);
        lon_sat2  = lon_sat1(kxy_indx);
        tim_sat2  = tim_sat1(kxy_indx);
        sen_sat2  = sen_sat1(kxy_indx);
        sol_sat2  = sol_sat1(kxy_indx);
        plv_sat2  = plv_sat1(:,kxy_indx);
        pwe_sat2  = pwe_sat1(:,kxy_indx);
        fgs_sat2  = fgs_sat1(:,kxy_indx);
        avk_sat2  = avk_sat1(:,kxy_indx);
%%% filter observation data based on their location (round 2 ~ 10km - round
%%% 3 ~ 1 km - round 1 ~ 100 km) or thinning process
        loc_r        = [lat_sat2',lon_sat2'];
        loc_rr       = round(loc_r,3);   
        [~,ia]       = unique(loc_rr,'rows','stable');
%%% =======================================================================
%%% if observation data is separately loaded         
%        lat_O       = nonzeros(lat_sat_all(:,step_size*(i-1)+nstep));
%        lon_O       = nonzeros(lon_sat_all(:,step_size*(i-1)+nstep));
%        loc_O       = [lat_O,lon_O];
%        loc_sat2    = [lat_sat2',lon_sat2'];
%        [~, idx_l]  = ismember(loc_O,loc_sat2, 'rows'); 
%        ia_l        = nonzeros(idx_l);
%        lat_sat  =  nonzeros(lat_sat_all(:,step_size*(i-1)+nstep));
%        lon_sat  =  nonzeros(lon_sat_all(:,step_size*(i-1)+nstep));
%        ch4_sat  =  ch4_sat_all(1:length(lat_sat),step_size*(i-1)+nstep);
%        sig_sat  =  sig_sat_all(1:length(lat_sat),step_size*(i-1)+nstep);
%        plv_sat  =  plv_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep);
%        pwe_sat  =  pwe_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep);
%        fgs_sat  =  fgs_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep);
%        avk_sat  =  avk_sat_all(:,1:length(lat_sat),step_size*(i-1)+nstep);
%      	 ch4_sat  = ch4_sat2(ia_l);
%        sig_sat  = sig_sat2(ia_l);
%        lat_sat  = lat_sat2(ia_l);
%        lon_sat  = lon_sat2(ia_l);
%        sen_sat  = sen_sat2(ia_l);
%        sol_sat  = sol_sat2(ia_l);
%        plv_sat  = plv_sat2(:,ia_l);
%        pwe_sat  = pwe_sat2(:,ia_l);
%        fgs_sat  = fgs_sat2(:,ia_l);
%        avk_sat  = avk_sat2(:,ia_l); 
%%% =======================================================================
        ch4_sat_bi  = ch4_sat2(ia);
        sig_sat  = sig_sat2(ia);
        lat_sat  = lat_sat2(ia);
        lon_sat  = lon_sat2(ia);
        tim_sat  = tim_sat2(ia);
        sen_sat  = sen_sat2(ia);
        sol_sat  = sol_sat2(ia);
        plv_sat  = plv_sat2(:,ia);
        pwe_sat  = pwe_sat2(:,ia);
        fgs_sat  = fgs_sat2(:,ia);
        avk_sat  = avk_sat2(:,ia); 
%%% bias correction (obsolete)        
       %ch4_sat  = ch4_sat_bi - (-0.46.*lat_sat - 13); % lattitude bias correction
       %ch4_sat  = ch4_sat_bi - (-0.55.*sol_sat - 12); % solar zenith angle bias correction
       %A_sat    = (1./cos(pi*sol_sat/180)) +  (1./cos(sen_sat*pi/180)); % AMF bais correction
       %ch4_sat  = ch4_sat_bi - (-21.*A_sat + 19); %AMF bais correction
%%% bias correction (latest update)
       %ch4_sat  = ch4_sat_bi - (-0.57.*lat_sat - 9.5); % lattitude bias correction
       %ch4_sat  = ch4_sat_bi - (-0.55.*sol_sat - 12); % solar zenith angle bias correction
       A_sat    = (1./cos(pi*sol_sat/180)) +  (1./cos(sen_sat*pi/180)); % AMF bais correction
       ch4_sat  = ch4_sat_bi - (-20.*A_sat + 17); %AMF bais correction
       %ch4_sat  = ch4_sat_bi;
   
   
%%% bias correction (latest update)
    	lat_sat_rad = deg2rad(lat_sat);
    	lon_sat_rad = deg2rad(lon_sat);
    	[x_sat,y_sat,z_sat] = sph2cart(lon_sat_rad,lat_sat_rad,r);
    	r_sat  = [x_sat(:),y_sat(:),z_sat(:)];
    	[xproj_sat,yproj_sat] = ll2psn(lat_sat,lon_sat,'TrueLat',45,'EarthRadius',6370000,...
    	'Eccentricity',1e-30,'meridian',-98);

%%% if no observations exist make variables empty  
    else
        disp('GOSAT Unavailable')
        lat_sat = [];
        lon_sat = [];
        tim_sat = [];
        sen_sat = [];
        sol_sat = [];
        ch4_sat = [];
        sig_sat = [];
        plv_sat = [];
        pwe_sat = [];
        fgs_sat = [];
        avk_sat = [];
    end
else
%%% if no observations exist make variables empty
disp('GOSAT Unavailable')
lat_sat = [];
lon_sat = [];
tim_sat = [];
sen_sat = [];
sol_sat = [];
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
