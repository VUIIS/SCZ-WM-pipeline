#!/usr/bin/env bash
#
# This script allows testing the compiled matlab, assuming the correct matlab
# runtime is installed. Much better to make sure it's working before actually 
# building the singularity container.

# Single fmri input
bin/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v97 \
    function scz_entrypoint \
    t1_niigz "/home/dylan/Documents/SCZ/INPUTS/t1.nii.gz" \
    fmri_niigz "{'/home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz'}" \
    out_dir "/home/dylan/Documents/SCZ/OUTPUTS" \
    xnat_project "TEST_PROJ" \
    xnat_subject "TEST_SUBJ" \
    xnat_session "TEST_SESS"

# Multiple fmri input
bin/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v97 \
    function scz_entrypoint \
    t1_niigz "/home/dylan/Documents/SCZ/INPUTS_v3/sub-OAS30001_ses-d0129_run-01_T1w.nii.gz" \
    fmri_niigz "{'/home/dylan/Documents/SCZ/INPUTS_v3/sub-OAS30001_ses-d0129_task-rest_run-02_bold.nii.gz', '/home/dylan/Documents/SCZ/INPUTS_v3/sub-OAS30001_ses-d0129_task-rest_run-03_bold.nii.gz'}" \
    out_dir "/home/dylan/Documents/SCZ/OUTPUTS" \
    xnat_project "TEST_PROJ" \
    xnat_subject "TEST_SUBJ" \
    xnat_session "TEST_SESS"