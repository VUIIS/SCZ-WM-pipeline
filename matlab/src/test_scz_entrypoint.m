% This script will test the matlab pipeline from the matlab command line -
% very useful for making sure it works before we bother to compile.

% SPM12 needs to be in the path also
addpath('/opt/spm12')
addpath(genpath('.'))

scz_entrypoint( ...
    't1_niigz','~/Documents/SCZ/INPUTS/t1.nii.gz', ...
	'fmri_niigz','~/Documents/SCZ/INPUTS/fmri.nii.gz', ...
	'out_dir','~/Documents/SCZ/OUTPUTS', ...
	'xnat_project','TEST_PROJ', ...
	'xnat_subject','TEST_SUBJ', ...
	'xnat_session','TEST_SESS', ...
	'xnat_scan_t1','TEST_SCAN_T1', ...
	'xnat_scan_fmri','TEST_SCAN_FMRI' ...
	);