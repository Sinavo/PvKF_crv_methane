%% 
%%% ======================================================================
%%% driver_hemi.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% ----------------------------------------------------------------------
%%%  NOTES:
%%%  This is the driver script for the PvKF methane data assimilation in 
%%%  hemispheric CMAQ domain using GOSAT observations. PvKF algorithm is 
%%%  based on a two-part paper (Voshtani et al., Remote Sensing, 2021 
%%%  , under review)
%%%  Other usefull materials:
%%%  doi:10.3402/tellusa.v68.31547
%%%  doi:10.3390/atmos9030086
%%%  doi:10.3390/atmos9020070
%%% =======================================================================

%%% Clear workspace and command window 
clc
clear all
close all

%%% Header
fprintf('\n ------------------------------------\n')
fprintf(' ---  Start of PvKF Assimilation  ---\n')
fprintf(' ------------------------------------\n')

%%% Read all input files
read_filename_hemi % check if this file exists in the current directory 
disp('succesfully read all inputs filenames')

%%% Some basic configurations
ncol = 187; 
nrow = 187;
nlev = 44;
B_step    = 24; % number of timestep in boundary file
step_size = 24; % number of timestep in a day  
end_day   = 1;  % last day of assimilation
start_day = 1;  % first day of assimilation

%%% Main loop of PvKF calculations 
main_loop_hemi_offline % check if this file exists in the current directory

%%% Finished simulation
fprintf('\n ------------------------------------\n')
fprintf(' ---    End of PvKF assimilation  ---\n')
fprintf(' ------------------------------------\n\n')
