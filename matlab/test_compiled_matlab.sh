#!/usr/bin/env bash
#
# This script allows testing the compiled matlab, assuming the correct matlab
# runtime is installed. Much better to make sure it's working before actually 
# building the singularity container.

bin/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v97 \
    t1_niigz /home/dylan/Documents/SCZ/INPUTS/t1.nii.gz \
    fmri_niigz /home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz \
    out_dir /home/dylan/Documents/SCZ/OUTPUTS