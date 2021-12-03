%% 
%%% ======================================================================
%%% main_loop_hemi_offline.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% PvKF for loop over days
for i = start_day:end_day
    tic
    
%%% extract daily values from netcdf files that are read 
    prep_day
    
%%% read observation data for the assimilation timestep    
    read_gosat
    %read_sciamachy
    
%%% to re-start from the middle
    if i == 1
       start_hour = 1;
     % load('RDA_w3f1q2_gos.mat')
       else
       start_hour = 1;
    end
    
%%% PvKF for loop over hours    
    for nstep=start_hour:1:step_size
   
%%% model forecast runs of concentration and error varaince        
%	run_hcmaq % check if this file exists in the current directory

%%% extract hourly values from netcdf files that are read 
	prep_hour_offline
    
%%% process satellite observation data: make it ready for assimilation     
    calc_gosat_offline
    %calc_sciamachy
    
%%% Analysis step of PvKF     
    analysis_offline

%%% replace analysis in model concentrations for the next forecast        
    write_cgrid_offline

save('gosat_cmaq_biascorrected_on_amf.mat','lat_sat_all','lon_sat_all','tim_sat_all','ch4_sat_all',...
     'sig_sat_all','sen_sat_all','sol_sat_all','obs_num_all','plv_sat_all','pwe_sat_all',...
     'fgs_sat_all','avk_sat_all','xm_o_avg_all')
    end
end
%%% =======================================================================
%%% END
%%% =======================================================================
