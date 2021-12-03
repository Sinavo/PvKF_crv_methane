%%
%%% ======================================================================
%%% write_cgrid.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

disp('writing analysis ...')

%%% replace analysis and analysis error variance for model forecast
%%% (nextstep)
%copyfile(con_name{step_size*(i-1)+nstep}, con_f_path, 'f');
%copyfile(var_name{step_size*(i-1)+nstep}, var_f_path, 'f');
netcdf.putVar(ncid_con,224,xa_3d);
netcdf.putVar(ncid_var,224,pa_3d);
netcdf.close(ncid_con)
netcdf.close(ncid_var)

%%% clear variables for next timestep
disp('cleaning memory ...')
clear  plv_sat pf_ii_m2o pres_m2o_l pwe_m pres_m2o_l indx_l hp_3d hcmaq_p
clear  xf_m2o pf_ii_m2o hp_m2o hp hph hp_o fgs_x pwe_sat
clear  xf_o_avg xa_hv xa_hv xf_m2o xf_o xf_o_avg pa_ii_hv
clear  avk_sat v_o v1 low_l pf_ii_hh pf_ii_hv hp_hv1

%%% =======================================================================
%%% END
%%% =======================================================================