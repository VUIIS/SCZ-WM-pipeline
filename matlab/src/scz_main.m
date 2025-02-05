% ----------------------added in v3 (start)----------------
% projID: project ID on XNAT (e.g., 'ADNI_23', 'BLSA', 'OASIS-3')
% ----------------------added in v3 (end)-----------------

% To add:
%   Support for multiple fmri/t1 files
%   Copy json with scan if it exists


function scz_main(inp)
t1_file = inp.t1_niigz;
fmri_file = inp.fmri_niigz;
xnat_project   = inp.xnat_project;
xnat_subject   = inp.xnat_subject;
xnat_session   = inp.xnat_session;
out_dir = inp.out_dir;

% Convert fmri_file to cell if multiple inputs are required
if contains(fmri_file,'}')
    eval(['fmri_file = ' fmri_file]);
else
    fmri_file=convertCharsToStrings(fmri_file);
end

% Convert chars to strings if deployed
if ischar(t1_file)
    t1_file=convertCharsToStrings(t1_file);
    xnat_project=convertCharsToStrings(xnat_project);
    xnat_subject=convertCharsToStrings(xnat_subject);
    xnat_session=convertCharsToStrings(xnat_session);
end

%% Prepare input files and move to output directory
disp('Preparing inputs')

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% preprocessing

% i=1;
% mkdir([out_dir '/preprocess/FunImg/' num2str(i)]);
% mkdir([out_dir '/preprocess/T1Img/' num2str(i)]);
% 
% tmp = split(fmri_file,'/');
% fmri_name = tmp{end};
% tmp = split(t1_file,'/');
% t1_name = tmp{end};
% %check for json
% tmp=dir(fmri_file);
% if isempty(dir([tmp.folder,'/*.json']))
%     fprintf('No .json files found. Extracting header info from nifti.');
% else
%     % copy json files
%     fmri_js = [fmri_file(1:end-7) '.json'];
%     copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),[fmri_name(1:end-7) '.json']));
%     t1_js = [t1_file(1:end-7) '.json'];
%     copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),[t1_name(1:end-7) '.json']));
% end
% copyfile(fmri_file,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),fmri_name)); % original T1 and rsfMRI
% copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i) ,t1_name));% original T1 and rsfMRI


if iscell(fmri_file)
    if contains(t1_file,'.gz')
        for i=1:length(fmri_file)
            
            mkdir([out_dir '/preprocess/FunImg/' num2str(i)]);
            mkdir([out_dir '/preprocess/T1Img/' num2str(i)]);
            
            
            tmp = split(fmri_file{i},'/');
            fmri_name = tmp{end};
            tmp = split(t1_file,'/');
            t1_name = tmp{end};
            %check for json
            tmp=dir(fmri_file{i});
            if isempty(dir([tmp.folder,'/*.json']))
                fprintf('No .json files found. Extracting header info from nifti.');
            else
                % copy json files
                fmri_js = [fmri_file{i}(1:end-7) '.json'];
                copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),[fmri_name(1:end-7) '.json']));
                t1_js = [t1_file(1:end-7) '.json'];
                copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),[t1_name(1:end-7) '.json']));
            end
            

            copyfile(fmri_file{i},sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),fmri_name)); % original T1 and rsfMRI
            copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),t1_name));% original T1 and rsfMRI
            
            gunzip(fullfile([out_dir '/preprocess/T1Img/' num2str(i)],t1_name));
            gunzip(fullfile([out_dir '/preprocess/FunImg/' num2str(i)],fmri_name));
        end
    else
        for i=1:length(fmri_file)
            mkdir([out_dir '/preprocess/FunImg/' num2str(i)]);
            mkdir([out_dir '/preprocess/T1Img/' num2str(i)]);
            
            
            tmp = split(fmri_file{i},'/');
            fmri_name = tmp{end};
            tmp = split(t1_file,'/');
            t1_name = tmp{end};
            %check for jsonn
            tmp=dir(fmri_file{i});
            if isempty(dir([tmp.folder,'/*.json']))
                fprintf('No .json files found. Extracting header info from nifti.');
            else
                % copy json files
                fmri_js = [fmri_file{i}(1:end-7) '.json'];
                copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),[fmri_name(1:end-4) '.json']));
                t1_js = [t1_file(1:end-7) '.json'];
                copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),[t1_name(1:end-4) '.json']));
            end
            
            fprintf(t1_file)
            fprintf(t1_name)
            fprintf(out_dir)
            
            copyfile(fmri_file{i},sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),fmri_name)); % original T1 and rsfMRI
            copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),t1_name));% original T1 and rsfMRI
        end
    end
else
    i=1;
    mkdir([out_dir '/preprocess/FunImg/' num2str(i)]);
    mkdir([out_dir '/preprocess/T1Img/' num2str(i)]);
    if contains(t1_file,'.gz')
        
        tmp = split(fmri_file,'/');
        fmri_name = tmp{end};
        tmp = split(t1_file,'/');
        t1_name = tmp{end};
        %check for json
        tmp=dir(fmri_file);
        if isempty(dir([tmp.folder,'/*.json']))
            fprintf('No .json files found. Extracting header info from nifti.');
        else
            % copy json files
            fmri_js = [fmri_file(1:end-7) '.json'];
            copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),[fmri_name(1:end-7) '.json']));
            t1_js = [t1_file(1:end-7) '.json'];
            copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),[t1_name(1:end-7) '.json']));
        end
        copyfile(fmri_file,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),fmri_name)); % original T1 and rsfMRI
        copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),t1_name));% original T1 and rsfMRI
        
        gunzip(fullfile([out_dir '/preprocess/T1Img/' num2str(i)],t1_name));
        gunzip(fullfile([out_dir '/preprocess/FunImg/' num2str(i)],fmri_name));
        
    else
        tmp = split(fmri_file,'/');
        fmri_name = tmp{end};
        tmp = split(t1_file,'/');
        t1_name = tmp{end};
        %check for json
        tmp=dir(fmri_file);
        if isempty(dir([tmp.folder,'/*.json']))
            fprintf('No .json files found. Extracting header info from nifti.');
        else
            % copy json files
            fmri_js = [fmri_file(1:end-7) '.json'];
            copyfile(fmri_js,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),[fmri_name(1:end-4) '.json']));
            t1_js = [t1_file(1:end-7) '.json'];
            copyfile(t1_js,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i),[t1_name(1:end-4) '.json']));
        end
        copyfile(fmri_file,sprintf('%s/preprocess/FunImg/%s/%s',out_dir,num2str(i),fmri_name)); % original T1 and rsfMRI
        copyfile(t1_file,sprintf('%s/preprocess/T1Img/%s/%s',out_dir,num2str(i) ,t1_name));% original T1 and rsfMRI
        
    end
end

%% ==== load config file =============================
% Change this path once code is containerized
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
    Cfg.SubjectID{1,1}=list(i).name;
    
    
    % get num of TRs and num of slices from fMRI header
    tmp = dir([out_dir,'/preprocess/FunImg/' list(i).name '/*.nii.gz']);
    info = niftiinfo([out_dir,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-3)]); % modified in v4: list(3) -> list(i) --------------------------------
    
    % ----------------------added in v3 (start)----------------
    pdim = info.PixelDimensions;
    if length(pdim)<4
        js = jsondecode(fileread([out_dir,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json'])); % modified in v4: list(3) -> list(i)
        TR = js.RepetitionTime;
    else
        TR = pdim(4);
    end
    disp(['TR =' num2str(TR) 's']);
    % ----------------------added in v3 (end)-----------------
    
    
    slnum = info.ImageSize(3); disp(['num of slices =' num2str(slnum)]);
    tpnum = info.ImageSize(4); disp(['num of time points =' num2str(tpnum)]);
    
    % ----------------------added in v3 (start)----------------
    % set default slice order according to the XNAT project
    if strcmp(xnat_project, 'ADNI_23')
        slorder = [1:2:slnum 2:2:slnum];
    elseif strcmp(xnat_project, 'BLSA')
        slorder = 1:slnum;
    elseif strcmp(xnat_project, 'OASIS-3')
        slorder = 1:slnum;
    else
        slorder = 1:slnum;
    end
    % check whether there is JSON file (BLSA-no; OASIS3-yes; ADNI23-most yes)
    % "SliceTiming" in JSON file indicates slice order
    % if there is no such info, then keep using default slice order
    if isfile([out_dir,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json']) % modified in v4: list(3) -> list(i) ---------------------------------------
        js = jsondecode(fileread([out_dir,'/preprocess/FunImg/' list(i).name '/' tmp.name(1:end-7) '.json'])); % modified in v4: list(3) -> list(i) ------------
        % check whether there is SliceTiming filed in JSON file (OASIS-part yes)
        if isfield(js, 'SliceTiming')
            slorder = round(js.SliceTiming/TR*slnum) + 1;
        end
    end
    % ----------------------added in v3 (end)-----------------
    
    %
    Cfg.TimePoints=tpnum; % set number of time points
    Cfg.TR=TR; % set TR
    Cfg.SliceTiming.SliceNumber=slnum; % set number of slices
    % ----------------------changed in v3 ---------------------
    Cfg.SliceTiming.SliceOrder=slorder; % set slice order % To Dylan: have to make sure
    % -------------------------------------------------------------
    
    Cfg.SliceTiming.ReferenceSlice=round(slnum/2); % set reference slice
    Cfg.IsCovremove=1; % do covariance removal
    Cfg.Covremove.CSF.IsRemove=1; % mean CSF as a regressor
    Cfg.Covremove.CSF.Mask='SPM'; % mask of CSF: SPM template
    
    
    
    
    %###############################################################%
    cd(out_dir);
    save([out_dir '/prepro_saved_' list(i).name '.mat'], 'Cfg', '-v7.3'); % modified in v4: prepro_saved.mat -> prepro_saved_subID.mat ---------------
    DPARSFA_run_XNAT([out_dir '/prepro_saved_' list(i).name '.mat']); % modified in v4: prepro_saved.mat -> prepro_saved_subID.mat; SCZ->XNAT ------------
    
end
cd(out_dir); 

% CHANGE TO OPTIONAL INPUT
wmthresh=0.80; % set white matter mask thresheold
gmthresh=0.80; % set gray  matter mask threshold

%% gen matrix
mkdir([out_dir '/result1_corrmatrix']);

copyfile([out_dir,'/preprocess/FunImgARCFWD/'],[out_dir,'/result1_corrmatrix/FunImgARCFWD/']);
copyfile([out_dir,'/preprocess/T1ImgNewSegment/'],[out_dir,'/result1_corrmatrix/T1ImgNewSegment/']);

%Get atlas path
tmp = which('Brodmann_YCG_newBA48.nii');
tmp=split(tmp,'/');
path = ['/', strjoin(tmp(2:end-3),'/')];

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
      
% --------------------- added in v2 -------------------------------------------------------------------------------------- 
% GM DKT ROI 
vv = spm_vol([path '/atlas/DKT/cropXm0Ym0Zm0_OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2.nii']); DKT = flip(spm_read_vols(vv)); 
DKT_label = unique(DKT(:)); 
DKT_label = DKT_label(2:end);
% GM HCPex ROI
vv = spm_vol([path '/atlas/HCPex/cropXm6Ym6Zm6_HCPex.nii']);  HCP = spm_read_vols(vv); 
HCP_label = 1:426; 
% GM Yeo-Schaefer-Kong 200 ROI 
vv = spm_vol([path '/atlas/Yeo_Schaefer_Kong/cropXm0Ym0Zm0_Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_1mm.nii']); SCH2 = flip(spm_read_vols(vv)); 
SCH2_label = 1:200;
% GM Yeo-Schaefer-Kong 400 ROI
vv = spm_vol([path '/atlas/Yeo_Schaefer_Kong/cropXm0Ym0Zm0_Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_1mm.nii']); SCH4 = flip(spm_read_vols(vv)); 
SCH4_label = 1:400;
%  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
       

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
    
        % ---------------------- added in v4 (start) ----------------------------------------------------------------------------------------------------------------------------
    % load configuration file of current subject (or acutally scan) 
    load([out_dir '/prepro_saved_' list(i).name '.mat']); 
    % ---------------------- added in v4 (end)  ----------------------------------------------------------------------------------------------------------------------------
                                
                            
    % init Brd_time / Eve_time / AAL time
    Brd_time = zeros(length(Brd_label),Cfg.TimePoints);
    Eve_time = zeros(length(Eve_label),Cfg.TimePoints);
    AAL_time = zeros(length(AAL_label),Cfg.TimePoints);
    
     % --------------------- added in v2 (start) ------------------------------------------------------------------------------------ 
    DKT_time = zeros(length(DKT_label),Cfg.TimePoints);
    HCP_time = zeros(length(HCP_label),Cfg.TimePoints);
    SCH2_time=zeros(length(SCH2_label),Cfg.TimePoints);
    SCH4_time=zeros(length(SCH4_label),Cfg.TimePoints);
    %  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
    
    
    % init corr matr
    matrEB = zeros(length(Eve_label), length(Brd_label));
    matrEA = zeros(length(Eve_label), length(AAL_label));
    
     % --------------------- added in v2 (start) ------------------------------------------------------------------------------------ 
    matrED = zeros(length(Eve_label), length(DKT_label));
    matrEH = zeros(length(Eve_label), length(HCP_label));   
    matrES2=zeros(length(Eve_label), length(SCH2_label));
    matrES4=zeros(length(Eve_label), length(SCH4_label));
    %  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
    
    
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
        roi_tc = img(logical(roi), :);
        Brd_time(j,:) = mean(roi_tc);
    end
    
    for j = 1:length(Eve_label)
        roi = logical(Eve==Eve_label(j)).*logical(wm>wmthresh);    
        roi_tc = img(logical(roi), :);
        Eve_time(j,:) = mean(roi_tc);
    end
    
    for j = 1:length(AAL_label)
        roi = logical(AAL==AAL_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(logical(roi), :);
        AAL_time(j,:) = mean(roi_tc); 
    end 
    
     % --------------------- added in v2 (start) ------------------------------------------------------------------------------------ 
    for j = 1:length(DKT_label)
        roi = logical(DKT==DKT_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(logical(roi), :);
        DKT_time(j,:) = mean(roi_tc); 
    end 
    for j = 1:length(HCP_label)
        roi = logical(HCP==HCP_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(logical(roi), :);
        HCP_time(j,:) = mean(roi_tc); 
    end 
    for j = 1:length(SCH2_label)
        roi = logical(SCH2==SCH2_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(logical(roi), :);
        SCH2_time(j,:) = mean(roi_tc); 
    end 
    for j = 1:length(SCH4_label)
        roi = logical(SCH4==SCH4_label(j)).*logical(gm>gmthresh);    
        roi_tc = img(logical(roi), :);
        SCH4_time(j,:) = mean(roi_tc); 
    end 
    %  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
    
    
    
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
    
    % --------------------- added in v2 (start) ------------------------------------------------------------------------------------
    for j = 1:length(Eve_label)
        for k = 1:length(DKT_label)
            tmp = corrcoef(Eve_time(j,:),DKT_time(k,:));
            matrED(j,k) = tmp(1,2);
        end
    end
    for j = 1:length(Eve_label)
        for k = 1:length(HCP_label)
            tmp = corrcoef(Eve_time(j,:),HCP_time(k,:));
            matrEH(j,k) = tmp(1,2);
        end
    end
    for j = 1:length(Eve_label)
        for k = 1:length(SCH2_label)
            tmp = corrcoef(Eve_time(j,:),SCH2_time(k,:));
            matrES2(j,k) = tmp(1,2);
        end
    end
    for j = 1:length(Eve_label)
        for k = 1:length(SCH4_label)
            tmp = corrcoef(Eve_time(j,:),SCH4_time(k,:));
            matrES4(j,k) = tmp(1,2);
        end
    end
    %  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
    
    
    save([out_dir,'/result1_corrmatrix/matr_', list(i).name, '.mat'], 'matrEB', 'matrEA', 'Brd_time','Eve_time', 'AAL_time');
    
    % --------------------- added in v2 (start) ------------------------------------------------------------------------------------ 
    save([out_dir,'/result1_corrmatrix/tc_', list(i).name, '.mat'], 'Brd_time','Eve_time', 'AAL_time', 'DKT_time', 'HCP_time', 'SCH2_time', 'SCH4_time');
    
    F = figure('visible','off');
    subplot(3,2,1); imagesc(matrEB);   axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-Brodmann');
    subplot(3,2,2); imagesc(matrEA);   axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-AAL');
    subplot(3,2,3); imagesc(matrED);   axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-DKT'); 
    subplot(3,2,4); imagesc(matrEH);   axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-HCPex'); 
    subplot(3,2,5); imagesc(matrES2); axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-Schaeler200'); 
    subplot(3,2,6); imagesc(matrES4); axis tight; colormap jet; caxis([-1 1]); colorbar; title('Eve-Schaeler400'); 
    
    saveas(gcf, [out_dir,'/result1_corrmatrix/matr_', list(i).name, '.png']);
    %  -------------------- added in v2 (end) --------------------------------------------------------------------------------------
    
    % --------------------- added by Dylan Lawless -----------------%
    % Make pdf
    mkdir([out_dir,'/PDF']);
    pdf_file = fullfile(out_dir,sprintf('/PDF/SCZ_WM_CorrMatr.pdf'));
    print(F,'-dpdf',pdf_file);
    
    % --------------------- added in v3 (start) ------------------------------------------------------------------------------------ 
    delete([out_dir '/result1_corrmatrix/FunImgARCFWD/',list(i).name,'/rDetrend_4DVolume.nii']);
    %  -------------------- added in v3 (end) --------------------------------------------------------------------------------------
    
end

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
    % ----------------------changed in v3 (start)--------------
    save([out_dir,'/result2_wm_alff/alff_' list(i).name(9:end-4) '.mat'],'alff');
    % ----------------------changed in v3 (end)---------------
    
    % --------------------- added in v3 (start) ------------------------------------------------------------------------------------ 
    delete([out_dir,'/result2_wm_alff/ALFF_FunImgARCFWD/r', list(i).name]);
    %  -------------------- added in v3 (end) --------------------------------------------------------------------------------------

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
    % ----------------------changed in v3 (start)--------------
    save([out_dir,'/result3_wm_falff/falff_' list(i).name(10:end-4) '.mat'],'falff');
    % ----------------------changed in v3 (end)---------------
    
    % --------------------- added in v3 (start) ------------------------------------------------------------------------------------ 
    delete([out_dir '/result3_wm_falff/fALFF_FunImgARCFWD/r' list(i).name]);
    %  -------------------- added in v3 (end) -------------------------------------------------------------------------------------

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
    % ----------------------changed in v3 (start)--------------
    save([out_dir '/result4_wm_reho/reho_' list(i).name(9:end-4) '.mat'],'reho');
    % ----------------------changed in v3 (end)---------------
    
    % --------------------- added in v3 (start) ------------------------------------------------------------------------------------ 
    delete([out_dir,'/result4_wm_reho/ReHo_FunImgARCFWD/r', list(i).name]);
    %  -------------------- added in v3 (end) -------------------------------------------------------------------------------------

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
%    save([out_dir,'/result5_wm_degree_centrality/degree_', xnat_subject, '.mat'],'degree');

    % ----------------------changed in v3 (start)--------------
    save([out_dir '/result5_wm_degree_centrality/degree_' list(i).name(46:end-4) '.mat'],'degree');
    % ----------------------changed in v3 (end)---------------
    
    % --------------------- added in v3 (start) ----------------------------------------------------------------------------------- 
    delete([out_dir,'/result5_wm_degree_centrality/DegreeCentrality_FunImgARCFWD/r', list(i).name]);
    %  -------------------- added in v3 (end) -------------------------------------------------------------------------------------

end

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

%% Zip all nifti files

cd(out_dir)

A = dir(out_dir);
for i = 3:length(A)
    if A(i).isdir == 0
        continue
    else 
        cd(A(i).folder);
        B = dir(A(i).name);
        for j = 3:length(B)
            if endsWith(B(j).name,'.nii')
                cd(B(j).folder);
                gzip(B(j).name);
                delete(B(j).name);
            elseif B(j).isdir == 0
                continue
            else
                cd(B(j).folder);
                C = dir(B(j).name);
                for k = 3:length(C)
                    if C(k).isdir == 1
                        cd(C(k).folder);
                        D = dir(C(k).name);
                        for l=3:length(D)
                            if D(l).isdir == 1
                                cd(D(l).folder);
                                E = dir(D(l).name);
                                for m=3:length(E)
                                    if endsWith(E(m).name,'.nii')
                                        cd(E(m).folder);
                                        gzip(E(m).name);
                                        delete(E(m).name);
                                    else
                                        continue
                                    end
                                end
                            elseif endsWith(D(l).name,'.nii')
                                cd(D(l).folder);
                                gzip(D(l).name);
                                delete(D(l).name);
                            else
                                continue
                            end
                        end
                    elseif endsWith(C(k).name,'.nii')
                        cd(C(k).folder);
                        gzip(C(k).name);
                        delete(C(k).name);
                    else
                        continue
                    end
                end
            end
        end
    end
end


disp('all done');
exit
end
