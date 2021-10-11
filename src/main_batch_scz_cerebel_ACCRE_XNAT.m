

function main_batch_scz_cerebel_ACCRE_XNAT(subIdx)

root = '/home/dylan/Documents/SCZ/';

addpath(genpath('/opt/spm12/'));
addpath(genpath([root 'z_batch_pipeline/add2path//NIFTI_20100819']));
addpath(genpath([root 'z_batch_pipeline/add2path/DPABI_V2.3_170105']));

list = dir([root 'data/']);
subID = list(subIdx+2).name;
disp(subID);

pwd = [root 'data/' subID]; % dir for current subject
path =[root 'z_batch_pipeline']; % dir for default cfg file


%% preprocessing
mkdir([pwd '/preprocess']); disp(pwd);

copyfile([pwd '/data_Funimg_T1img/*'],[pwd '/preprocess/']); % original T1 and rsfMRI

load([path '/prepro_saved.mat']); % load default cfg file
%######################## set parameters #######################%
% set dir etc.
Cfg.WorkingDir=[pwd,'/preprocess'];
Cfg.DataProcessDir=[pwd,'/preprocess'];
list=dir([pwd,'/preprocess/FunImg/']);
Cfg.SubjectID={};
for i=3:length(list)
    Cfg.SubjectID{i-2,1}=list(i).name;
end

% set TR (nii.hdr.dime.pixdim(5) is not always TR)
TR = 2; disp(['TR =' num2str(TR) 's']); % To Dylan: might be different case by case
% get num of TRs and num of slices from fMRI header
tmp = dir([pwd,'/preprocess/FunImg/' list(3).name '/*.nii.gz']);
nii = load_untouch_nii_gz([pwd,'/preprocess/FunImg/' list(3).name '/' tmp.name]); 
slnum = nii.hdr.dime.dim(4); disp(['num of slices =' num2str(slnum)]);
tpnum = nii.hdr.dime.dim(5); disp(['num of time points =' num2str(tpnum)]);
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

wmthresh=0.80; % set white matter mask thresheold
gmthresh=0.80; % set gray  matter mask threshold
%###############################################################%
save([pwd '/prepro_saved.mat'], 'Cfg', '-v7.3'); % save to subject's dir
DPARSFA_run_SCZ([pwd '/prepro_saved.mat']); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cd(pwd);  % load([pwd '/prepro_saved.mat']);

%% gen matrix
mkdir([pwd '/result1_corrmatrix']);

copyfile([pwd,'/preprocess/FunImgARCFWD/'],[pwd,'/result1_corrmatrix/FunImgARCFWD/']);
copyfile([pwd,'/preprocess/T1ImgNewSegment/'],[pwd,'/result1_corrmatrix/T1ImgNewSegment/']);
% GM Brodmann ROI
nii=load_nii([path,'/atlas/Brodmann/Brodmann_YCG_newBA48.nii']); Brd=nii.img;
load([path,'/atlas/Brodmann/Brodmann_YCG_Labels.mat']); label=cell2mat(Reference(:,2)); Brd_label=label([2:2:82,83:-2:3]);
% WM Eve ROI
nii=load_nii([path,'/atlas/EVE+Cerebellar/Eve+Cerebellar.nii']); Eve=nii.img;
Eve_label=[ 8 10 12  2 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 ...
            3  4  5  6 ...
           47 45 43 41 39 37 35 33 31 29 27 25 23 21 19 17 15 13  1 11 9 7];
% GM AAL ROI
nii=load_nii([path '/atlas/AAL/AAL3v1_1mm_improv.nii']); AAL=nii.img;
AAL_label = 1:170;
       
list=dir([pwd,'/result1_corrmatrix/FunImgARCFWD/']); list(1:2)=[];
for i=1:length(list)
    % reslice Fun image (61x73x61x300 -> 181x217x181)
    reslice_nii([pwd,'/result1_corrmatrix/FunImgARCFWD/',list(i).name,'/Detrend_4DVolume.nii'],...
                [pwd,'/result1_corrmatrix/FunImgARCFWD/',list(i).name,'/image.nii'],[1 1 1]);  
    % reslice GM image (121x145x121 -> 181x217x181)
    tmp=dir([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/wc1*.nii']); % gray matter segment
    reslice_nii([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/',tmp(1).name],...
                [pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/image_noflip_gm.nii'],[1 1 1]);
    % reslice WM image (121x145x121 -> 181x217x181) 
    tmp=dir([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/wc2*.nii']); % white matter segment
    reslice_nii([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/',tmp(1).name],...
                [pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/image_noflip_wm.nii'],[1 1 1]);
    
    % init Brd_time / Eve_time / AAL time
    Brd_time=zeros(length(Brd_label),Cfg.TimePoints);
    Eve_time=zeros(length(Eve_label),Cfg.TimePoints);
    AAL_time=zeros(length(AAL_label),Cfg.TimePoints);
    
    % init corr matr
    matrEB=zeros(length(Eve_label), length(Brd_label));
    matrEA=zeros(length(Eve_label), length(AAL_label));
    
    % load resliced fun / gm / wm
    nii=load_nii([pwd,'/result1_corrmatrix/FunImgARCFWD/',list(i).name,'/image.nii']);      img=nii.img;
    nii=load_nii([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/image_noflip_gm.nii']);gm =nii.img;
    nii=load_nii([pwd,'/result1_corrmatrix/T1ImgNewSegment/',list(i).name,'/image_noflip_wm.nii']);wm =nii.img; % % 
    
    img=reshape(img,181*217*181,Cfg.TimePoints);
    imgstd=repmat(std(img,[],2),1,Cfg.TimePoints); % %
    imgmean=repmat(mean(img,2),1,Cfg.TimePoints);
    img=(img-imgmean)./(imgstd+eps);
    img=reshape(img,181,217,181,Cfg.TimePoints);

    for j=1:length(Brd_label)
        roi=logical(Brd==Brd_label(j)).*logical(gm>gmthresh);
        roi=repmat(roi,1,1,1,Cfg.TimePoints);
        tmp=img.*roi;
        tmp=reshape(tmp,181*217*181,Cfg.TimePoints);
        roi=reshape(roi,181*217*181,Cfg.TimePoints);
        Brd_time(j,:)=sum(tmp,1)./sum(roi,1);
    end
    
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);    
        roi=repmat(roi,1,1,1,Cfg.TimePoints);
        tmp=img.*roi;
        tmp=reshape(tmp,181*217*181,Cfg.TimePoints);
        roi=reshape(roi,181*217*181,Cfg.TimePoints);
        Eve_time(j,:)=sum(tmp,1)./sum(roi,1);
    end
    
    for j=1:length(AAL_label)
        roi=logical(AAL==AAL_label(j)).*logical(gm>gmthresh);    
        roi=repmat(roi,1,1,1,Cfg.TimePoints);
        tmp=img.*roi;
        tmp=reshape(tmp,181*217*181,Cfg.TimePoints);
        roi=reshape(roi,181*217*181,Cfg.TimePoints);
        AAL_time(j,:)=sum(tmp,1)./sum(roi,1);
     end 
    
    for j=1:length(Eve_label)
        for k=1:length(Brd_label)
            tmp=corrcoef(Eve_time(j,:),Brd_time(k,:));
            matrEB(j,k)=tmp(1,2);
        end
    end

    for j=1:length(Eve_label)
        for k=1:length(AAL_label)
            tmp=corrcoef(Eve_time(j,:),AAL_time(k,:));
            matrEA(j,k)=tmp(1,2);
        end
    end
%     matr=mapminmax(matr(:)');
%     matr=reshape(matr,length(Eve_label),length(Brd_label));
    save([pwd,'/result1_corrmatrix/matr2_',list(i).name,'_noflip.mat'],'matrEB', 'matrEA', 'Brd_time','Eve_time', 'AAL_time');
end
list1=list;

%% gen wm alff
mkdir([pwd '/result2_wm_alff']);

copyfile([pwd,'/preprocess/Results/ALFF_FunImgARCFWD'],[pwd,'/result2_wm_alff/ALFF_FunImgARCFWD']);
list=dir([pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/ALFF*.nii']);
for i=1:length(list)
    name=[pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/',list(i).name];
    tmpname=list(i).name(8:end);
    newname=[pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/image',tmpname];
    reslice_nii(name,newname,[1 1 1]);
    nii=load_nii([pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/image',tmpname]);
    img=nii.img;      img=img(end:-1:1,:,:,:);
    nii.img=img;
    save_nii(nii,[pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/image',tmpname]);
    nii=load_nii([pwd,'/result2_wm_alff/ALFF_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    alff=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        alff(j)=sum(tmp(:))./sum(roi(:));
    end
    tmpname2=tmpname(1:end-4);
    save([pwd,'/result2_wm_alff/alff',tmpname2,'.mat'],'alff');
end

%% gen wm falff
mkdir([pwd '/result3_wm_falff']);

copyfile([pwd,'/preprocess/Results/fALFF_FunImgARCFWD'],[pwd,'/result3_wm_falff/fALFF_FunImgARCFWD']);
list=dir([pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/fALFF*.nii']);
for i=1:length(list)
    name=[pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/',list(i).name];
    tmpname=list(i).name(9:end);
    newname=[pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/image',tmpname];
    reslice_nii(name,newname,[1 1 1]);
    nii=load_nii([pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    img=img(end:-1:1,:,:,:);
    nii.img=img;
    save_nii(nii,[pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/image',tmpname]);
    nii=load_nii([pwd,'/result3_wm_falff/fALFF_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    falff=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        falff(j)=sum(tmp(:))./sum(roi(:));
    end
    tmpname2=tmpname(1:end-4);
    save([pwd,'/result3_wm_falff/falff',tmpname2,'.mat'],'falff');
end

%% gen wm reho
mkdir([pwd '/result4_wm_reho']);

copyfile([pwd,'/preprocess/Results/ReHo_FunImgARCFWD'],[pwd,'/result4_wm_reho/ReHo_FunImgARCFWD']);
list=dir([pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/ReHo*.nii']);
for i=1:length(list)
    name=[pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/',list(i).name];
    tmpname=list(i).name(8:end);
    newname=[pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/image',tmpname];
    reslice_nii(name,newname,[1 1 1]);
    nii=load_nii([pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    img=img(end:-1:1,:,:,:);
    nii.img=img;
    save_nii(nii,[pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/image',tmpname]);
    nii=load_nii([pwd,'/result4_wm_reho/ReHo_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    reho=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        reho(j)=sum(tmp(:))./sum(roi(:));
    end
    tmpname2=tmpname(1:end-4);
    save([pwd,'/result4_wm_reho/reho',tmpname2,'.mat'],'reho');
end

%% gen wm degree_centrality
mkdir([pwd '/result5_wm_degree_centrality']);

copyfile([pwd,'/preprocess/Results/DegreeCentrality_FunImgARCFWD'],[pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD']);
list=dir([pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/DegreeCentrality_PositiveWeighted*.nii']);
for i=1:length(list)
    name=[pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/',list(i).name];
    tmpname=list(i).name(45:end);
    newname=[pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/image',tmpname];
    reslice_nii(name,newname,[1 1 1]);
    nii=load_nii([pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    img=img(end:-1:1,:,:,:);
    nii.img=img;
    save_nii(nii,[pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/image',tmpname]);
    nii=load_nii([pwd,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/image',tmpname]);
    img=nii.img;
    degree=zeros(1,length(Eve_label));
    for j=1:length(Eve_label)
        roi=logical(Eve==Eve_label(j)).*logical(wm>wmthresh);
        tmp=img.*roi;
        degree(j)=sum(tmp(:))./sum(roi(:));
    end
    tmpname2=tmpname(1:end-4);
    save([pwd,'/result5_wm_degree_centrality/degree',tmpname2,'.mat'],'degree');
end


%% copy some important files
copyfile([pwd,'/preprocess/PicturesForChkNormalization/*'],[pwd,'/']);
copyfile([pwd,'/preprocess/RealignParameter/HeadMotion.mat'],[pwd,'/']);


%% remove preprocess to save some space in ACCRE
if exist([pwd '/preprocess/T1Img/'], 'dir'), rmdir([pwd '/preprocess/T1Img/'], 's'); end
if exist([pwd '/preprocess/T1ImgCoreg/'], 'dir'), rmdir([pwd '/preprocess/T1ImgCoreg/'], 's'); end
if exist([pwd '/preprocess/T1ImgNewSegment/'], 'dir'), rmdir([pwd '/preprocess/T1ImgNewSegment/'], 's'); end
if exist([pwd '/preprocess/FunImg/'], 'dir'), rmdir([pwd '/preprocess/FunImg/'], 's'); end
if exist([pwd '/preprocess/FunImgA/'], 'dir'), rmdir([pwd '/preprocess/FunImgA/'], 's'); end
if exist([pwd '/preprocess/FunImgAR/'], 'dir'), rmdir([pwd '/preprocess/FunImgAR/'], 's'); end
if exist([pwd '/preprocess/FunImgARC/'], 'dir'), rmdir([pwd '/preprocess/FunImgARC/'], 's'); end
if exist([pwd '/preprocess/FunImgARCF/'], 'dir'), rmdir([pwd '/preprocess/FunImgARCF/'], 's'); end
if exist([pwd '/preprocess/FunImgARCFW/'], 'dir'), rmdir([pwd '/preprocess/FunImgARCFW/'], 's'); end
if exist([pwd '/preprocess/FunImgARCFWD/'], 'dir'), rmdir([pwd '/preprocess/FunImgARCFWD/'], 's'); end
if exist([pwd '/result1_corrmatrix/FunImgARCFWD/',list1(i).name,'/image.nii'], 'file') 
    delete([pwd '/result1_corrmatrix/FunImgARCFWD/',list1(i).name,'/image.nii']);
end

disp('all done');
% exit
end
