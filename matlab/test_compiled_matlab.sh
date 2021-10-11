#!/usr/bin/env bash
#
# This script allows testing the compiled matlab, assuming the correct matlab
# runtime is installed. Much better to make sure it's working before actually 
# building the singularity container.

bin/run_matlab_entrypoint.sh /usr/local/MATLAB/MATLAB_Runtime/v97 \
    t1_niigz ../INPUTS/t1.nii.gz \
    seg_niigz ../INPUTS/seg.nii.gz \
    diameter_mm 30 \
    out_dir ../OUTPUTS