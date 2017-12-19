%% demo code for voxel-wise HRF deconvolution without finer temporal grid, but with the option of choosing the estimation type (FLAG)
%% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
%% Guo-Rong Wu, gronwu@gmail.com, UESTC, UGent, 2013.9.12
%% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D..
%% A blind deconvolution approach to recover effective connectivity brain networks
%% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .
clc,clear
matlabpool open 'local2' 6
addpath(pwd)
% add path to SPM here if not already there
%% open Matlab parallel computing, NumWorkers: set a reasonable number yourself. If you don't have parallel facilities no prob, but change "parfor" to normal "for"
% try
%     myCluster = parcluster('local');
%     myCluster.NumWorkers = 8;
%     saveAsProfile(myCluster,'local2');
%     matlabpool open 'local2' 8
% end
TR = .72; %in seconds
thr = 1; % threshold, for example 1 SD.
event_lag_max_seconds=9;    % the (estimated) maximum lagged time from neural event to BOLD event, in seconds.
event_lag_max =round(9/TR); % the (estimated) maximum lagged time from neural event to BOLD event, in points.
T = round(30/TR); %this is the maximum length of the whole lifespan of the HRF, assumed to be 30 seconds
maskfile='mask_3mm_278ROIs.nii'; %path to the mask with the cortical voxels
main= '/home/daniele';  % main root directory, change here if you don't run the script from the directory where it's located.
%data_dir = fullfile(main,'FunImgNormalizedCovremovedDetrendedFiltered'); %% where data are stored after smoothing,regression, filtering, detrending or whatever preprocessing
data_dir=fullfile(main,'/RS_HCP'); %directory containing the folders of each subject, sub1, sub2 etc
dir_out=main; %change also here. This is the directory to store the results, best use the previous directory, going a folder back
sub = dir(data_dir);
sub(1:2)=[]; %this removes the "." and ".."
flag_HRF=3; % flag 1 for FIR, 2 for canonical, 3 for rbeta
dir_hrf={'/FIR/','/canon/','/rbeta/'};
save_dir = fullfile(dir_out,'/results_deconv_HCP',dir_hrf{flag_HRF}); %% save dir, change the name as you like.

% make directories for the HRF parameters
mkdir(save_dir);

save pars % here you save the initialization parameters

for isub=1%:length(sub)
    disp('Reading data ...')
    sub_dir = fullfile(data_dir,sub(isub).name);
    disp(sub(isub).name)
    cd(sub_dir);
    clear imag
    % if your preprocessed data are not stored in an image, but in a
    % vector, you can call this vector rsig and skip the following lines
    imag = dir('*.nii'); %% if *.nii, change it yourself.
    tic
    
    if length(imag)>1
        fourD=false;
        [data1] = spm_read_vols(spm_vol(maskfile));
        data1(isnan(data1))=0;
        voxel_ind=find(data1);
        num_voxel=length(voxel_ind);
        rsig = zeros(size(imag,1),num_voxel);
        
        parfor k = 1:length(imag)
            [data1] = spm_read_vols(spm_vol(imag(k).name));
            rsig(k,:) =  data1(voxel_ind);
        end
        rsig=single(rsig);clear data1
    else
        fourD=true;
        data1=spm_read_vols(spm_vol(imag.name));
        data1(isnan(data1))=0;
        nimag=size(data1,4);
        data1=reshape(data1,[],nimag);
        i1 = spm_read_vols(spm_vol(maskfile));
        voxel_ind=find(i1);
        num_voxel=length(voxel_ind);
        datamask=data1(voxel_ind,:);
        rsig=datamask';clear datamask
        rsig=single(rsig);
    end
    toc
    disp('Done,')
    disp('Filtering and Detrending ...')
        tic
        rsig = rest_IdealFilter(rsig, TR, [0.01 0.1]); % bandpass filter
        rsig = spm_detrend(rsig,3); % ensure stability
        toc
        disp('Done')
    disp('Retrieving HRF ...');
    tic
    %[data_deconv onset hrf event_lag PARA] = wgr_deconv_canonhrf_par(rsig,thr,event_lag_max,TR);
    [sig_deconv, hdrf] = hrf_retrieval_and_deconvolution_para(rsig,thr,event_lag_max,TR,T,flag_HRF);
    toc
    disp('Done');
    tic
    disp('Saving ..');
    save(fullfile(save_dir,[sub(isub).name,'_mask_hrf.mat']),'hdrf','sig_deconv','-v7.3');
    toc
    %%
    % if you want to write the deconvolved data as images, use
    % the following lines, otherwise comment
    tic
    disp('writing back deconvolved data to .nii images ...')
    v=spm_vol(maskfile);
    v.dt=[16,0];
    sub_save_dir = fullfile(save_dir,sub(isub).name);
    mkdir(sub_save_dir)
    % writing back into nifti files
    if fourD
        for k = 1:nimag
            v.fname = fullfile(sub_save_dir,['imag',num2str(k,'%04d'),'.nii']);
            data = zeros(v.dim);
            data(voxel_ind) = sig_deconv(k,:);
            spm_write_vol(v,data);
        end
    else
        for k = 1:length(imag)
            v.fname = fullfile(sub_save_dir,imag(k).name);
            data = zeros(v.dim);
            data(voxel_ind) = sig_deconv(k,:);
            spm_write_vol(v,data);
        end
    end
    toc
    disp('Done');
    
    %%
    
    
% here in order to free some memory you clear everything and reload the
% parameters that you saved before entering the loop
    clear 
     load('/media/daniele/driveD/hrf/pars');
end