#!/usr/bin/env bash
#
# Main pipeline
echo Running $(basename "${BASH_SOURCE}")

# Execute matlab command
run_spm12.sh "${MATLAB_RUNTIME}" function scz_entrypoint \
    t1_niigz "${t1_niigz}" \
    fmri_niigz "${fmri_niigz}" \
    out_dir "${out_dir}" \
    xnat_project "${xnat_project}" \
    xnat_subject "${xnat_subject}" \
    xnat_session "${xnat_session}"

# PDF creation and zipping nifti files
cd ${out_dir}
make_pdf.py

