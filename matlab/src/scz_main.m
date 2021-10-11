% is the config file necessary?
% change way niftis are read

function scz_main(inp)
t1_file = inp.t1_niigz;
fmri_file = inp.fmri_niigz;
xnat_project   = inp.xnat_project;
xnat_subject   = inp.xnat_subject;
xnat_session   = inp.xnat_session;
xnat_scan_t1   = inp.xnat_scan_t1;
xnat_scan_fmri = inp.xnat_scan_fmri;
out_dir = inp.out_dir;

fprintf('%s %s %s - T1 %s, fMRI %s\n', ...
	xnat_project,xnat_subject,xnat_session, ...
	xnat_scan_t1,xnat_scan_fmri);
fprintf('t1_niigz:   %s\n', t1_file   );
fprintf('fmri_niigz: %s\n', fmri_file );
fprintf('out_dir:    %s\n', out_dir   );


% Nifti read/write is handled this way:
% https://github.com/VUIIS/spm_readwrite_nii

% do I need to do this or get the full SPM?

%% Prepare input files and move to output directory
disp('Preparing inputs')

%use out_path instead of root

% Add to path: SPM12, DPABI
%addpath(genpath('/Users/dylanlawless/Documents/spm12'));
addpath(genpath('DPABI_V2.3_170105'));

%pwd = [root 'data/' subID]; % dir for current subject
path = '/opt/scz/bin';
addpath(genpath('/opt/scz/bin'));


%% preprocessing
mkdir([out_dir '/preprocess/FunImg/' xnat_subject]);
mkdir([out_dir '/preprocess/T1Img/' xnat_subject]);

copyfile(fmri_file,fullfile([out_dir '/preprocess/FunImg/' xnat_subject],'fmri.nii.gz')); % original T1 and rsfMRI
copyfile(t1_file,fullfile([out_dir '/preprocess/T1Img/' xnat_subject],'t1.nii.gz'));% original T1 and rsfMRI

gunzip(fullfile([out_dir '/preprocess/T1Img/' xnat_subject],'t1.nii.gz'));
gunzip(fullfile([out_dir '/preprocess/FunImg/' xnat_subject],'fmri.nii.gz'));

fmri_nii = fullfile([out_dir '/preprocess/FunImg/' xnat_subject],'fmri.nii');
t1_nii = fullfile([out_dir '/preprocess/T1Img/' xnat_subject],'t1.nii');

%% ==== load config file =============================

LOADFILENAME1=which(fullfile('prepro_saved.mat'));

% Load the data file from current working directory
load(which(LOADFILENAME1),'Cfg');

%% ######################## set parameters #######################%
% set dir etc.
Cfg.WorkingDir=[out_dir,'/preprocess'];
Cfg.DataProcessDir=[out_dir,'/preprocess'];
list=dir([out_dir,'/preprocess/FunImg/']);
Cfg.SubjectID={};
for i=3:length(list)
    Cfg.SubjectID{i-2,1}=list(i).name;
end

% get num of TRs and num of slices from fMRI header

info = niftiinfo(fmri_nii);
TR = info.PixelDimensions(4); disp(['TR =' num2str(TR) 's']);
slnum = info.ImageSize(3); disp(['num of slices =' num2str(slnum)]);
tpnum = info.ImageSize(4); disp(['num of time points =' num2str(tpnum)]);
clear nii;
%
Cfg.TimePoints=tpnum; % set number of time points
Cfg.TR=TR; % set TR
Cfg.SliceTiming.SliceNumber=slnum; % set number of slices
Cfg.SliceTiming.SliceOrder=[1:2:slnum 2:2:slnum]; % set slice order % To Dylan: have to make sure
Cfg.SliceTiming.ReferenceSlice=round(slnum/2); % set reference slice
Cfg.IsCovremove=1; % do covariance removal
Cfg.Covremove.CSF.IsRemove=1; % mean CSF as a regressor
Cfg.Covremove.CSF.Mask='SPM'; % mask of CSF: SPM template

% CHANGE TO OPTIONAL INPUT
wmthresh=0.80; % set white matter mask thresheold
gmthresh=0.80; % set gray  matter mask threshold


%###############################################################%
cd(out_dir); 
save([out_dir '/prepro_saved.mat'], 'Cfg', '-v7.3'); % save to subject's dir
DPARSFA_run_SCZ([out_dir '/prepro_saved.mat']); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cd(out_dir); 

%% gen matrix
mkdir([out_dir '/result1_corrmatrix']);

copyfile([out_dir,'/preprocess/FunImgARCFWD/'],[out_dir,'/result1_corrmatrix/FunImgARCFWD/']);
copyfile([out_dir,'/preprocess/T1ImgNewSegment/'],[out_dir,'/result1_corrmatrix/T1ImgNewSegment/']);
% GM Brodmann ROI
vv = spm_vol([path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii']); Brd = spm_read_vols(vv);
load([path,'/atlas/Brodmann/Brodmann_YCG_Labels.mat']); label=cell2mat(Reference(:,2)); Brd_label=label([2:2:82,83:-2:3]);
% WM Eve ROI
vv = spm_vol([path,'/atlas/EVE+Cerebellar/Eve+Cerebellar.nii']); Eve = spm_read_vols(vv);
Eve_label = [ 8 10 12  2 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 ...
                    3  4  5  6 ...
                    47 45 43 41 39 37 35 33 31 29 27 25 23 21 19 17 15 13  1 11 9 7];
% GM AAL ROI
vv = spm_vol([path '/atlas/AAL/AAL3v1_1mm_improv.nii']); AAL = spm_read_vols(vv);
AAL_label = 1:170;
       
list = dir([out_dir,'/result1_corrmatrix/FunImgARCFWD/']); list(1:2)=[];
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
for i = 1:length(list)
    % reslice Fun image (61x73x61x300 -> 181x217x181)
    spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result1_corrmatrix/FunImgARCFWD/',list(i).name,'/Detrend_4DVolume.nii']}, flags);
                
    % reslice GM image (121x145x121 -> 181x217x181)
    tmp = dir([out_dir,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/wc1*.nii']); % gray matter segment
    spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/',tmp(1).name]}, flags);              
                
    % reslice WM image (121x145x121 -> 181x217x181) 
    tmp = dir([out_dir,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/wc2*.nii']); % white matter segment
    spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/',tmp(1).name]}, flags);
    
    % init Brd_time / Eve_time / AAL time
    Brd_time = zeros(length(Brd_label),Cfg.TimePoints);
    Eve_time = zeros(length(Eve_label),Cfg.TimePoints);
    AAL_time = zeros(length(AAL_label),Cfg.TimePoints);
    
    % init corr matr
    matrEB = zeros(length(Eve_label), length(Brd_label));
    matrEA = zeros(length(Eve_label), length(AAL_label));
    
    % load resliced fun / gm / wm
    vv = spm_vol([out_dir,'/result1_corrmatrix/FunImgARCFWD/', list(i).name, '/rDetrend_4DVolume.nii']); img = spm_read_vols(vv);
    tmp = dir([out_dir,'/result1_corrmatrix/T1ImgNewSegment/', list(i).name,'/rwc1*.nii']);
    vv = spm_vol([out_dir,'/result1_corrmatrix/T1ImgNewSegment/', list(i).name, '/', tmp(1).name]); gm = spm_read_vols(vv);
    tmp = dir([out_dir,'/result1_corrmatrix/T1ImgNewSegment/', list(i).name,'/rwc2*.nii']);
    vv = spm_vol([out_dir,'/result1_corrmatrix/T1ImgNewSegment/', list(i).name, '/', tmp(1).name]); wm = spm_read_vols(vv); 
    % normalize time courses   
    img = reshape(img,181*217*181,Cfg.TimePoints); % N x nTimePoint
    nrows = 181*217*181;
    for ii = 1:nrows
        tc = img(ii, :); 
        if sum(abs(tc))==0, continue, end
        img(ii, :) = (tc - mean(tc))/(std(tc) + eps);  
    end

    for j = 1:length(Brd_label)
        roi = logical(Brd==Brd_label(j)).*logical(gm>gmthresh);
        roi_tc = img(roi, :);
        Brd_time(j,:) = mean(roi_tc);
    end
    
    for j = 1:length(Eve_label)
        roi = logical(Eve==Eve_label(j)).*logical(wm>wmthresh);    
        roi_tc = img(roi, :);
        Eve_time(j,:) = mean(roi_tc);
    end
    
    for j = 1:length(AAL_label)
        roi = logical(AAL==AAL_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(roi, :);
        AAL_time(j,:) = mean(roi_tc); 
     end 
    
    for j = 1:length(Eve_label)
        for k = 1:length(Brd_label)
            tmp = corrcoef(Eve_time(j,:),Brd_time(k,:));
            matrEB(j,k) = tmp(1,2);
        end
    end

    for j = 1:length(Eve_label)
        for k = 1:length(AAL_label)
            tmp = corrcoef(Eve_time(j,:),AAL_time(k,:));
            matrEA(j,k) = tmp(1,2);
        end
    end
    
    save([out_dir,'/result1_corrmatrix/matr_', list(i).name, '.mat'], 'matrEB', 'matrEA', 'Brd_time','Eve_time', 'AAL_time');
end
list1=list;

%% gen wm alff
mkdir([out_dir '/result2_wm_alff']);

copyfile([out_dir,'/preprocess/Results/ALFF_FunImgARCFWD'], [out_dir,'/result2_wm_alff/ALFF_FunImgARCFWD']);
list = dir([out_dir,'/result2_wm_alff/ALFF_FunImgARCFWD/ALFF*.nii']);
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
for i = 1:length(list)  
    spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result2_wm_alff/ALFF_FunImgARCFWD/', list(i).name]}, flags);      
    vv = spm_vol([out_dir,'/result2_wm_alff/ALFF_FunImgARCFWD/r', list(i).name]); img = spm_read_vols(vv);

    alff = zeros(1,length(Eve_label));
    for j = 1:length(Eve_label)
        roi = logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp = img.*roi;
        alff(j) = sum(tmp(:))./sum(roi(:));
    end
    save([out_dir,'/result2_wm_alff/alff_', subID, '.mat'],'alff');
end

%% gen wm falff
mkdir([out_dir '/result3_wm_falff']);

copyfile([out_dir,'/preprocess/Results/fALFF_FunImgARCFWD'],[out_dir,'/result3_wm_falff/fALFF_FunImgARCFWD']);
list=dir([out_dir,'/result3_wm_falff/fALFF_FunImgARCFWD/fALFF*.nii']);
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
for i=1:length(list)
    spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result3_wm_falff/fALFF_FunImgARCFWD/', list(i).name]}, flags);
    vv = spm_vol([out_dir,'/result3_wm_falff/fALFF_FunImgARCFWD/r', list(i).name]); img = spm_read_vols(vv);    

    falff=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        falff(j)=sum(tmp(:))./sum(roi(:));
    end
    save([out_dir,'/result3_wm_falff/falff_', subID, '.mat'],'falff');
end

%% gen wm reho
mkdir([out_dir '/result4_wm_reho']);

copyfile([out_dir,'/preprocess/Results/ReHo_FunImgARCFWD'],[out_dir,'/result4_wm_reho/ReHo_FunImgARCFWD']);
list=dir([out_dir,'/result4_wm_reho/ReHo_FunImgARCFWD/ReHo*.nii']);
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
for i=1:length(list)
     spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result4_wm_reho/ReHo_FunImgARCFWD/', list(i).name]}, flags);  
    vv = spm_vol([out_dir,'/result4_wm_reho/ReHo_FunImgARCFWD/r', list(i).name]); img = spm_read_vols(vv); 

    reho=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        reho(j)=sum(tmp(:))./sum(roi(:));
    end
    save([out_dir,'/result4_wm_reho/reho_', subID, '.mat'],'reho');
end

%% gen wm degree_centrality
mkdir([out_dir '/result5_wm_degree_centrality']);

copyfile([out_dir,'/preprocess/Results/DegreeCentrality_FunImgARCFWD'],[out_dir,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD']);
list=dir([out_dir,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/DegreeCentrality_PositiveWeighted*.nii']);
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
for i=1:length(list)
     spm_reslice_quiet({[path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii'] ...
                                [out_dir,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/', list(i).name]}, flags);  
    vv = spm_vol([out_dir,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/r', list(i).name]); img = spm_read_vols(vv);  

    degree=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        degree(j)=sum(tmp(:))./sum(roi(:));
    end
    save([out_dir,'/result5_wm_degree_centrality/degree_', subID, '.mat'],'degree');
end


%% copy some important files
copyfile([out_dir,'/preprocess/PicturesForChkNormalization/*'],[out_dir,'/']);
copyfile([out_dir,'/preprocess/RealignParameter/HeadMotion.mat'],[out_dir,'/']);


%% remove preprocess to save some space in ACCRE
if exist([out_dir '/preprocess/T1Img/'], 'dir'), rmdir([out_dir '/preprocess/T1Img/'], 's'); end
if exist([out_dir '/preprocess/T1ImgCoreg/'], 'dir'), rmdir([out_dir '/preprocess/T1ImgCoreg/'], 's'); end
if exist([out_dir '/preprocess/T1ImgNewSegment/'], 'dir'), rmdir([out_dir '/preprocess/T1ImgNewSegment/'], 's'); end
if exist([out_dir '/preprocess/FunImg/'], 'dir'), rmdir([out_dir '/preprocess/FunImg/'], 's'); end
if exist([out_dir '/preprocess/FunImgA/'], 'dir'), rmdir([out_dir '/preprocess/FunImgA/'], 's'); end
if exist([out_dir '/preprocess/FunImgAR/'], 'dir'), rmdir([out_dir '/preprocess/FunImgAR/'], 's'); end
if exist([out_dir '/preprocess/FunImgARC/'], 'dir'), rmdir([out_dir '/preprocess/FunImgARC/'], 's'); end
if exist([out_dir '/preprocess/FunImgARCF/'], 'dir'), rmdir([out_dir '/preprocess/FunImgARCF/'], 's'); end
if exist([out_dir '/preprocess/FunImgARCFW/'], 'dir'), rmdir([out_dir '/preprocess/FunImgARCFW/'], 's'); end
if exist([out_dir '/preprocess/FunImgARCFWD/'], 'dir'), rmdir([out_dir '/preprocess/FunImgARCFWD/'], 's'); end
if exist([out_dir '/result1_corrmatrix/FunImgARCFWD/',list1(i).name,'/rDetrend_4DVolume.nii'], 'file') 
    delete([out_dir '/result1_corrmatrix/FunImgARCFWD/',list1(i).name,'/rDetrend_4DVolume.nii']);
end

disp('all done');
exit
end
