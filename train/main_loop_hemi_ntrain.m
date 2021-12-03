%% 
%%% ======================================================================
%%% main_loop_hemi_ntrain.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% comment "load" below if bias-corrected GOSAT data is not used
load('cv_input_train23_all.mat','lat_sat_all_train'...
     ,'lon_sat_all_train','tim_sat_all_train','ch4_sat_all_train',...
     'sig_sat_all_train','plv_sat_all_train','pwe_sat_all_train',...
     'fgs_sat_all_train','avk_sat_all_train')
%%% cv_input_train12_all.mat cv_input_train13_all.mat (3-fold
%%% cross-validation require 3 input of bias corrected train data

%%% PvKF for loop over days
for i = start_day:end_day
    tic
    
    %%% extract daily values from netcdf files that are read 
    prep_day
    %read_gosat
    %read_sciamachy

%%% to re-start from the middle
    if i == 1
       start_hour = 1;
     %  load('train_Lh350_Lv1_i2f45q20.mat')
       else
       start_hour = 1;
    end
    
    for nstep=start_hour:1:step_size

%%% model forecast runs of concentration and error varaince        
%	run_hcmaq % check if this file exists in the current directory
	
%%% extract hourly values from netcdf files that are read 
	prep_hour % check if this file exists in the current directory
    
%%% process satellite observation data: make it ready for assimilation              
        calc_gosat_crossv          
        %calc_sciamachy

%%% Analysis step of PvKF        
    analysis_plus

%%% replace analysis in model concentrations for the next forecast    
    write_cgrid
    
%%% save important variables for the post-processing
save('ntrain_Lh400_Lv2_i2f45q20.mat','lat_sat_all','lon_sat_all','ch4_sat_all','sig_sat_all',...
     'plv_sat_all','pwe_sat_all','fgs_sat_all','avk_sat_all',...
     'chi_sq_all','obs_num_all','innov_all','xf_o_avg_all','var_B_all')
    end
end

%%% =======================================================================
%%% END
%%% =======================================================================
