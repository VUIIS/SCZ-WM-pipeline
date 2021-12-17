#!/usr/bin/env bash
#
# Test the pipeline outside the container. Be sure the src directory is in the
# path.

export MATLAB_RUNTIME=/usr/local/MATLAB/MATLAB_Runtime/v97


export t1_niigz="/home/dylan/Documents/SCZ/INPUTS/t1.nii.gz"
export fmri_niigz="{'/home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz'}"
export xnat_project="TEST_PROJ"
export xnat_session="TEST_SESS"
export xnat_subject="TEST_SUBJ"
export out_dir="/home/dylan/Documents/SCZ/OUTPUTS"

pipeline_entrypoint.sh \
    --t1_niigz "${ti_niigz}" \
    --fmri_niigz "${fmri_niigz}" \
    --out_dir "${out_dir}" \
    --xnat_project "${xnat_project}" \
    --xnat_subject "${xnat_subject}" \
    --xnat_session "${xnat_session}"
