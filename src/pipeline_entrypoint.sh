#!/usr/bin/env bash
#
# Primary entrypoint for our pipeline. This just parses the command line 
# arguments, exporting them in environment variables for easy access
# by other shell scripts later. Then it calls the rest of the pipeline.
#
# Example usage:
# 
# pipeline_entrypoint.sh --image_niigz /path/to/image.nii.gz --diameter_mm 30

# This statement at the top of every bash script is helpful for debugging
echo Running $(basename "${BASH_SOURCE}")

# Initialize defaults for any input parameters where that seems useful
export label_info="UNKNOWN SCAN"
export out_dir=/OUTPUTS

# Parse input options
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        
        --t1_niigz)
            # Our example code takes a single 3D T1 nifti image as input. This
            # is expected to be the fully qualified path and filename.
            export t1_niigz="$2"; shift; shift ;;

        --fmri_niigz)
            # Segmentation of the T1. Useful here to provide example code for
            # viewing ROIs. Expected to be in the same geometry, position, voxel
            # size, etc - e.g. output of slant or multi-atlas pipelines.
            export seg_niigz="$2"; shift; shift ;;

        --label_info)
            # Labels from XNAT that we will use to label the QA PDF
            export label_info="$2"; shift; shift ;;

        --out_dir)
            # Where outputs will be stored. Also the working directory
            export out_dir="$2"; shift; shift ;;

        *)
            echo "Input ${1} not recognized"
            shift ;;

    esac
done


# Now that we have all the inputs stored in environment variables, call the
# main pipeline. We run it in xvfb so that we have a virtual display available.
xvfb-run -n $(($$ + 99)) -s '-screen 0 1600x1200x24 -ac +extension GLX' \
    bash pipeline_main.sh
