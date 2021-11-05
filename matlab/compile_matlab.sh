#!/bin/bash

#
# Compile the matlab code so we can run it without a matlab license. To create a 
# linux container, we need to compile on a linux machine. That means a VM, if we 
# are working on OS X.
#
# We require on our compilation machine:
#     Matlab 2019b, including compiler, with license
#     Installation of SPM12, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
#
# The matlab version matters. If we compile with R2019b, it will only run under 
# the R2019b Runtime.
#
# The SPM12 version also matters. The compilation code is written for r7771.


# Where to find SPM12 on our compilation machine
SPM_PATH=/opt/spm12

# We may need to add Matlab to the path on the compilation machine
export PATH=/usr/local/MATLAB/R2019b/bin:${PATH}

# We use SPM12's standalone tool, but edited to add our own code to the 
# compilation path
WD=`pwd`
matlab -nodisplay -nodesktop -nosplash -sd "${WD}" -r \
    "spm_make_standalone_local('${SPM_PATH}','${WD}/bin','${WD}/src'); exit"

# We grant lenient execute permissions to the matlab executable and runscript so
# we don't have hiccups later.
chmod go+rx "${WD}"/bin/spm12
chmod go+rx "${WD}"/bin/run_spm12.sh