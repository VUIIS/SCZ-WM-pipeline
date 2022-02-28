#!/usr/bin/env bash
#
# Test the built singularity container.

export t1_niigz="/home/dylan/Documents/SCZ/INPUTS/t1.nii.gz"
export fmri_niigz="/home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz"
export xnat_project="TEST_PROJ"
export xnat_session="TEST_SESS"
export xnat_subject="TEST_SUBJ"
export out_dir="/home/dylan/Documents/SCZ/OUTPUTS"


singularity run --cleanenv --contain \
    --home $(pwd -P) \
    --bind $(pwd -P)/INPUTS:/INPUTS \
    --bind $(pwd -P)/OUTPUTS:/OUTPUTS \
    scz_wm_v2.simg \
    --t1_niigz "${t1_niigz}" \
    --fmri_niigz "${fmri_niigz}" \
    --out_dir "${out_dir}" \
    --xnat_project "${xnat_project}" \
    --xnat_subject "${xnat_subject}" \
    --xnat_session "${xnat_session}"

