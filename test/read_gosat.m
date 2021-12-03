%% 
%%% ======================================================================
%%% read_gosat.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% Read GOSAT data variables
lat_draw   = double(ncread(Oname_go{i},'latitude'));
lon_draw   = double(ncread(Oname_go{i},'longitude'));
tim_draw   = double(ncread(Oname_go{i},'time'));
ch4_draw   = double(ncread(Oname_go{i},'xch4'));
sig_draw   = double(ncread(Oname_go{i},'xch4_uncertainty'));
avk_draw   = double(ncread(Oname_go{i},'xch4_averaging_kernel'));
plv_draw   = double(ncread(Oname_go{i},'pressure_levels'));
pwe_draw   = double(ncread(Oname_go{i},'pressure_weight'));
fgs_draw   = double(ncread(Oname_go{i},'ch4_profile_apriori'));
sen_draw   = double(ncread(Oname_go{i},'sensor_zenith_angle'));
sol_draw   = double(ncread(Oname_go{i},'solar_zenith_angle'));

%%% Read flags and create index based on (change for user flags) 
idx_quality  = logical(ncread(Oname_go{i},'xch4_quality_flag')); % good quality flag
idx_landtype = logical(ncread(Oname_go{i},'flag_landtype')); % land flag
idx_sunglint = logical(ncread(Oname_go{i},'flag_sunglint'));% sunglint flag
idx_ql       = bsxfun(@and,~idx_quality,~idx_landtype);

%%% select indexed variables 
lat_ql   = lat_draw(idx_ql);
lon_ql   = lon_draw(idx_ql);
sen_ql   = sen_draw(idx_ql);
sol_ql   = sol_draw(idx_ql);
tim_ql   = tim_draw(idx_ql);
ch4_ql   = ch4_draw(idx_ql);
sig_ql   = sig_draw(idx_ql);
avk_ql   = avk_draw(:,idx_ql);
plv_ql   = plv_draw(:,idx_ql);
pwe_ql   = pwe_draw(:,idx_ql);
fgs_ql   = fgs_draw(:,idx_ql);
indx_nan   = isnan(ch4_ql);
lat        = lat_ql(~indx_nan);
lon        = lon_ql(~indx_nan);
sen        = sen_ql(~indx_nan);
sol        = sol_ql(~indx_nan);
tim        = tim_ql(~indx_nan);
avk        = avk_ql(:,~indx_nan);
plv        = plv_ql(:,~indx_nan);
pwe        = pwe_ql(:,~indx_nan);
fgs        = fgs_ql(:,~indx_nan); fgs_d(isnan(fgs_ql)) = 1800;
ch4_d      = ch4_ql(~indx_nan);
sig_d      = sig_ql(~indx_nan); sig_d = fillmissing(sig_ql,'previous');
ch4 = filloutliers(ch4_d,'nearest','mean');
sig = filloutliers(sig_d,'nearest','mean');

%%% convert pressure hPa to Pa
plv   = plv *100;

%%% time conversion for the delivered data 
t        = datetime(tim,'ConvertFrom','posixtime');
t1       = datetime(2010,04,i+1,0,0,0);
t_step   = dateshift(t1,'start','hour',0:24/step_size:24);
t_step2   = dateshift(t1,'start','hour',0:24/B_step:24);
hm       = 100*hour(t) + minute(t);
hm_start = 100*hour(t1) + minute(t1);
hm_step  = 100*hour(t_step) + minute(t_step);
hm_step2  = 100*hour(t_step2) + minute(t_step2);

%%% =======================================================================
%%% END
%%% =======================================================================