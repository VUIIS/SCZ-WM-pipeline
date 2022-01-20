# SCZ White Matter Pipeline

This spider will preprocess fMRI data as well as corresponding T1 data, extract mean time-courses of each predefined ROI and compute the correlation matrices between white matter ROIs and gray matter ROIs. Please see Gao’s publications [1, 2] for more details. The spider will also compute FALFF, ALFF and ReHo maps. 

This XNAT spider is currently designed for three databases (ADNI_23, BLSA and OASIS-3) which are proposed to be analyzed in white matter reanalysis project (PI: Dr. Gore and Dr. Landman). 

## Inputs: 
fMRI (.nii.gz)

T1 (.nii.gz)

Configuration file (.mat)
  
## Output:
Preprocessed fMRI in MNI space: 

	../scz_OUTPUTS/result1_corrmatrix/FunImgARCFWD/1/Detrend_4DVolume.nii.gz

Tissue probability maps (gray matter and white matter) in MNI space:

 	../scz_OUTPUTS/result1_corrmatrix/T1ImgNewSegment/1/wc1T1.nii.gz

	../scz_OUTPUTS/result1_corrmatrix/T1ImgNewSegment/1/wc2T1.nii.gz

Functional connectivity matrices between white matter ROIs and gray matter ROIs: 

	../scz_OUTPUTS/result1_corrmatrix/matr_1.mat

	../scz_OUTPUTS/result1_corrmatrix/matr_1.png

Mean time-courses of the white and gray matter ROIs:

	../ scz_OUTPUTS/result1_corrmatrix/tc_1.mat	

Whole-brain ALFF maps:

	../ scz_OUTPUTS/result2_wm_alff/ALFF_FunImgARCFWD/ALFFMap_1.nii.gz

	../ scz_OUTPUTS/result2_wm_alff/ALFF_FunImgARCFWD/mALFFMap_1.nii.gz

	../ scz_OUTPUTS/result2_wm_alff/ALFF_FunImgARCFWD/zALFFMap_1.nii.gz

Whole-brain FALFF maps:

	../ scz_OUTPUTS/result3_wm_falff/fALFF_FunImgARCFWD/fALFFMap_1.nii.gz

	../ scz_OUTPUTS/result3_wm_falff/fALFF_FunImgARCFWD/mfALFFMap_1.nii.gz

	../ scz_OUTPUTS/result3_wm_falff/fALFF_FunImgARCFWD/zfALFFMap_1.nii.gz

Whole brain ReHo maps:

	../ scz_OUTPUTS/result4_wm_reho/ReHo_FunImgARCFWD/ReHoMap_1.nii.gz

	../ scz_OUTPUTS/result4_wm_reho/ReHo_FunImgARCFWD/mReHoMap_1.nii.gz

	../ scz_OUTPUTS/result4_wm_reho/ReHo_FunImgARCFWD/zReHoMap_1.nii.gz


Whole brain maps of degree of centrality:

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/DegreeCentrality_PositiveBinarizedSumBrainMap_1.nii.gz

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/DegreeCentrality_PositiveWeightedSumBrainMap_1.nii.gz

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/mDegreeCentrality_PositiveBinarizedSumBrainMap_1.nii.gz

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/mDegreeCentrality_PositiveWeightedSumBrainMap_1.nii.gz

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/zDegreeCentrality_PositiveBinarizedSumBrainMap_1.nii.gz

	../ scz_OUTPUTS/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/zDegreeCentrality_PositiveWeightedSumBrainMap_1.nii.gz

# References
[1] Gao Y, Sengupta A, Li M, et al. (2020) Functional connectivity of white matter as a biomarker of cognitive decline in Alzheimer’s disease. PLoS ONE 15(10): e0240513. https://doi.org/10.1371/journal.pone.0240513

[2] Gao Y, Li M, Huang AS. Lower functional connectivity of white matter during rest and working memory tasks is associated with cognitive impairments in schizophrenia. Schizophr Res. 2021 Jul;233:101-110. doi: 10.1016/j.schres.2021.06.013. Epub 2021 Jun 29. PMID: 34215467; PMCID: PMC8442250. 

