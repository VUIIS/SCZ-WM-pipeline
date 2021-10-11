#!/usr/bin/env bash
#
# Test the built singularity container.


singularity run --cleanenv --contain \
    --home $(pwd -P)/INPUTS \
    --bind workdir:/tmp \
    --bind INPUTS:/INPUTS \
    --bind OUTPUTS:/OUTPUTS \
    demo.simg \
    --t1_niigz /INPUTS/t1.nii.gz \
    --fmri_niigz /INPUTS/seg.nii.gz \
    --diameter_mm 30 \
    --out_dir /OUTPUTS

