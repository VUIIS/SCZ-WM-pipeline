function [meanfmri_niigz,rfmri_niigz,srfmri_niigz] = zip_outputs( ...
	meanfmri_nii,rfmri_nii,srfmri_nii)

% Zip images
gzip(meanfmri_nii);
meanfmri_niigz = [meanfmri_nii '.gz'];

gzip(rfmri_nii);
rfmri_niigz = [rfmri_nii '.gz'];

gzip(srfmri_nii);
srfmri_niigz = [srfmri_nii '.gz'];