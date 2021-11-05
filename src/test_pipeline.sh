#!/usr/bin/env bash
#
# Test the pipeline outside the container. Be sure the src directory is in the
# path.

# Just the PDF creation part
export t1_niigz="/home/dylan/Documents/SCZ/INPUTS/t1.nii.gz"
export fmri_niigz="/home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz"
export xnat_project="TEST_PROJ"
export xnat_session="TEST_SESS"
export xnat_subject="TEST_SUBJ"
export out_dir="/home/dylan/Documents/SCZ/OUTPUTS"


# The entire thing
pipeline_entrypoint.sh \
    --t1_niigz /home/dylan/Documents/SCZ/INPUTS/t1.nii.gz \
    --fmri_niigz /home/dylan/Documents/SCZ/INPUTS/fmri.nii.gz \
    --out_dir /home/dylan/Documents/SCZ/OUTPUTS
