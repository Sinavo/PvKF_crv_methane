%% 
%%% ======================================================================
%%% read_filename_hemi.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% read the names of concentrations netcdf output files
con_path = '/home/sinavo/pkg/cmaq/5.3/data/cctm_hemi/hemi_h1/con9';
con_file = dir(fullfile(con_path,'CCTM_CGRID_*'));
con_name = arrayfun(@(X) fullfile(con_path, X.name), con_file, 'uni', 0);
disp('reading concentration filenames ...')

%%% read the names of variances netcdf output files
var_path = '/home/sinavo/pkg/cmaq/5.3/data/cctm_hemi/hemi_h1/var9';
var_file = dir(fullfile(var_path,'CCTM_CGRID_*'));
var_name = arrayfun(@(X) fullfile(var_path, X.name), var_file, 'uni', 0);
disp('reading varaince filenames ...')

%%% read the names of meteorology netcdf (MCIP) output files
Mpath          = '/home/sinavo/pkg/cmaq/5.3/data/mcip_hemi/mcip_h_5d';
Mfile_cro2d    = dir(fullfile(Mpath,'METCRO2D*'));
Mfile_cro3d    = dir(fullfile(Mpath,'METCRO3D*'));
Mname_cro2d    = arrayfun(@(X) fullfile(Mpath, X.name), Mfile_cro2d...
    , 'uni', 0);
Mname_cro3d    = arrayfun(@(X) fullfile(Mpath, X.name), Mfile_cro3d...
    , 'uni', 0);

%%% read the names of lat-lon netcdf (MCIP) output files
Lpath          = '/home/sinavo/pkg/cmaq/5.3/data/mcip_hemi/mcip_h_5d';
Lfile          = dir(fullfile(Lpath,'GRIDCRO2D*'));
Lname          = arrayfun(@(X) fullfile(Lpath, X.name), Lfile, 'uni', 0);
disp('reading met-data filename ...')

%%% read the names of observations (GOSAT, SCIAMACHY, IASI,...) files
Opath_go       = '/home/sinavo/pkg/cmaq/5.3/data/obs/gosat';
% Opath_sc     = '/home/sinavo/pkg/cmaq/5.3/data/obs/sciamachy';
% Opath_ia     = '/home/sinavo/pkg/cmaq/5.0.2/data_sina/obs/iasi';
% Opath_ar     = '/home/sinavo/pkg/cmaq/5.0.2/data_sina/obs/airs';
% Opath_te     = '/home/sinavo/pkg/cmaq/5.0.2/data_sina/obs/tes';
Ofile_go       = dir(fullfile(Opath_go,'ESACCI-GHG-L2-CH4-GOSAT-SRPR*'));
% Ofile_sc     = dir(fullfile(Opath_sc,'ESACCI-GHG-L2-CH4-SCIAMACHY-I*'));
% Ofile_ia     = dir(fullfile(Opath_ia,'CH4_IASIA_NLIS*'));
% Ofile_te     = dir(fullfile(Opath_te,'TES-Aura_2010.04*.he5'));
% Ofile_ar_std = dir(fullfile(Opath_ar,'AIRS*Std*hdf'));
% Ofile_ar_sup = dir(fullfile(Opath_ar,'AIRS*SUB*hdf'));
Oname_go       = arrayfun(@(X) fullfile(Opath_go, X.name), Ofile_go...
    , 'uni', 0);
% Oname_sc     = arrayfun(@(X) fullfile(Opath_sc, X.name), Ofile_sc...
%, 'uni', 0);
%Oname_ia     = arrayfun(@(X) fullfile(Opath_ia, X.name), Ofile_ia...
%, 'uni', 0);
%Oname_te     = arrayfun(@(X) fullfile(Opath_te, X.name), Ofile_te...
%, 'uni', 0);
%Oname_ar_std = arrayfun(@(X) fullfile(Opath_ar, X.name), Ofile_ar_std...
%, 'uni', 0);
%Oname_ar_sup = arrayfun(@(X) fullfile(Opath_ar, X.name), Ofile_ar_sup...
%, 'uni', 0);
disp('reading observation filenames ...')

%%% =======================================================================
%%% END
%%% =======================================================================
