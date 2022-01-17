# PvKF assimilation (+Cross-Validation)
Readme file for parametric variance Kalman filter (PvKF: https://github.com/Sinavo/PvKF_reg_methane.git) with cross-validation capability from Voshtani et al., (2022), Part II, DOI: https://doi.org/10.3390/rs14020375 <br /> 
PvKF is demonstrated for GOSAT methane assimilation with hemispheric CMAQ model <br />
December 1, 2021 <br />

# info

- Cross-validation is based on 3-fold approach. For GOSAT observations sample outputs of "train or active" observations (cv_input_train23_all.mat) and "test or passive" observations (cv_input_test1_all.mat) are generated. Users can also generate their own specific k-fold cross-validation train and test sets by running the "cv_make_3fold.m" file. <br />
- In the train folder, you can run PvKF assimilation with the training set (2/3 of observations data as a default). The output (.mat) and concentrations/variances NetCDF files are then used to run the PvKF in the test folder. Running PvKF assimilation in both train and test folders is the same as regular PvKF: https://github.com/Sinavo/PvKF_reg_methane.git.<br /> 
- The assimilation can be run using the "driver_hemi.m" file in Matlab or using the "run_driver_hemi.sh" file in bash in both train and test folders.<br /> 
- The assimilation uses publicly available datasets and models, including <br />

CMAQ (https://doi.org/10.5281/zenodo.1212601)<br />
SMOKE (https://doi.org/10.5281/zenodo.4480334)<br />
GOSAT CH4 Proxy (https://catalogue.ceda.ac.uk/uuid/96d5b75ea29946c5aab8214ddbab252b)<br />
GOSAT CH4 Full Physics (https://catalogue.ceda.ac.uk/uuid/f9154243fd8744bdaf2a59c39033e659)<br />
EDGAR CH4 (https://edgar.jrc.ec.europa.eu/index.php/dataset_ghg60)

NOTE: GOSTA bias-corrected data against surface (NOAA) observations are generated as "gosat_cmaq_biascorrected_qc_on_latitude.mat". Bias correction is explained in Voshtani et al., (2022, Part I, Section 4.1 & 4.2). The code will still run without loading the bias-corrected (.mat) file, once it is commented. Input and modelling data can also be accessed by contacting Sina Voshtani (sinavoshtani@cmail.carleton.ca)




