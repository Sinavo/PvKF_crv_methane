#!/bin/bash
module load matlab
matlab -nodesktop -nojvm  -r "run /home/sinavo/pkg/cmaq/5.3/CCTM/da_9/mat_hemi/driver_hemi.m;"
