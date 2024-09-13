#!/usr/bin/env bash
#
# Primary entrypoint for our pipeline. This just parses the command line 
# arguments, exporting them in environment variables for easy access
# by other shell scripts later. Then it calls the rest of the pipeline.
#
# Example usage:
# 
# pipeline_entrypoint.sh --t1_niigz image.nii.gz
#                        --fmri_niigz {image1.nii.gz,image2.nii.gz}
#                        --out_dir /output/directory
#                        --xnat_project XNAT_project_name
#                        --xnat_subject XNAT_subject_name
#                        --xnat_session XNAT_session_name
#

echo Running $(basename "${BASH_SOURCE}")

# Initialize defaults 
export xnat_subject="TEST_SUBJ"
export out_dir=/OUTPUTS

# Parse input options
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        
        --t1_niigz)
            export t1_niigz="$2"; shift; shift ;;

        --fmri_niigz)
            export fmri_niigz="$2"; shift; shift ;;

        --out_dir)
            export out_dir="$2"; shift; shift ;;

        --xnat_project)
            export xnat_project="$2"; shift; shift ;;

        --xnat_subject)
            export xnat_subject="$2"; shift; shift ;;

        --xnat_session)
            export xnat_session="$2"; shift; shift ;;

        *)
            echo "Input ${1} not recognized. Exiting program."
            exit
            shift ;;

    esac
done

# Run pipeline with xvfb for display window
xvfb-run -n $(($$ + 99)) -s '-screen 0 1600x1200x24 -ac +extension GLX' \
    bash pipeline_main.sh
