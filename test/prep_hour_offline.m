%% 
%%% ======================================================================
%%% prep_hour_offline.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% read concentration netcdf variable
ncid_con      = netcdf.open(con_name{step_size*(i-1)+nstep},'WRITE');
c_con         = double(netcdf.getVar(ncid_con,224));
%c_con         = double(netcdf.getVar(ncid_con,1));

%%% read variance netcdf variable
ncid_var      = netcdf.open(var_name{step_size*(i-1)+nstep},'WRITE');
ncid_var_0    = netcdf.open(var_name{1},'WRITE');

%%% modelling error factor (fq is in Part I, Eq.41)
%c_var_0       = double(netcdf.getVar(ncid_var_0,224));
c_var_0       = double(netcdf.getVar(ncid_var_0,1));
fq            = 1;
%qq             = (c_var_0 .^ 0.5) ./ (fq*24);  

%%% add a uniform (additive) model error to the forecast error variance
%c_var_t       = double(netcdf.getVar(ncid_var,224));
c_var_t       = double(netcdf.getVar(ncid_var,1));
%c_var         = c_var_t + qq.^2;  
c_var         = c_var_t;  

%%% compute layer pressures based on the surface pressur and sigma value 
hcmaq_p_s    = squeeze(hcmaq_p_s_25(:,:,:,nstep));
ptop_hcmaq   = 5000* ones(size(hcmaq_p_s));
for l=1:45
    hcmaq_p(:,:,l)  = ( hcmaq_siglvl(l) * (hcmaq_p_s - ptop_hcmaq) ) + ptop_hcmaq ;
end
%%% =======================================================================
%%% END
%%% =======================================================================
