%% demo code for voxel-wise HRF deconvolution
%% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
%% Guo-Rong Wu, gronwu@gmail.com, UESTC, UGent, 2013.9.12
%% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D..
%% A blind deconvolution approach to recover effective connectivity brain networks
%% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .
clc,clear;warning off all

%% Mask file
%%===================================
brainmask = 'mask_3mm_278ROIs.nii';
%%===================================
v=spm_vol(brainmask);
brain = spm_read_vols(v);
voxel_ind = find(brain==11); %% change as you like. brain>0 for all brain, or brain==ROIid
num_voxel = length(voxel_ind);
v.dt=[16,0];

%% open Matlab parallel computing, NumWorkers: set a reasonable number yourself. If you don't have parallel facilities no prob, but change "parfor" to normal "for"
try
    myCluster = parcluster('local');
    myCluster.NumWorkers = 8;
    saveAsProfile(myCluster,'local2');
    parpool open 'local2' 8
end

%%===========PARAMETERS========================
% choose the set of basis functions THIS MUST BE AN INPUT

%para.estimation='canon2dd'; % this one for canonical HRF plus two derivatives
para.estimation='sFIR'; % this one for smoothed FIR
%para.estimation='FIR'; % this one for unsmoothed FIR

temporal_mask = []; % without mask, it means temporal_mask = ones(nobs,1); i.e. all time points included. nobs: number of observation = size(data,1). if want to exclude the first 1~5 time points, let temporal_mask(1:5)=0;

TR = 1.5; % THIS WILL BE READ FROM THE BIDS DATA

para.TR = TR;

%%% the following parameter (upsample grid) can be > 1 only for Canonical.
%%% Set = 1 for FIR
para.T  = 3; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid

para.T0 = 3; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
if para.T==1
    para.T0 = 1;
end

min_onset_search = 4; % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = 8; % maximum delay allowed between event and HRF onset (seconds)

para.dt  = para.TR/para.T; % fine scale time resolution.

para.TD_DD = 2; % time and dispersion derivative

para.AR_lag = 1; % AR(1) noise autocorrelation.

para.thr = 1; % (mean+) para.thr*standard deviation threshold to detect event.

para.len = 24; % length of HRF, in seconds

para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);

%%===================================

%%===========fMRI Data========================
main= pwd;  % data directory
save_dir = fullfile(main,'output_HRF_Deconv'); %% save dir, change the name as you like.

%% list all your 4D files (change names and directories accordingly)
sub={
    [pwd '\input_BOLD_data\sub_001\processed_fmri_data_4D.nii']
    };

for isub=1:length(sub)
    disp('Reading data ...')
    [sub_dir,name,ext] = fileparts(sub{isub,1});
    sub_save_dir = strrep(sub_dir,main,save_dir);
    mkdir(sub_save_dir)
    
    v1 = spm_vol(sub{isub,1});
    if ~all(v1(1).dim==v.dim)
        error('The dimension of your mask is different than the one of your fMRI data!')
    end
    data1 = spm_read_vols(v1);
    nobs = size(data1,4);
    data1 = reshape(data1,[],nobs)';
    bold_sig =  data1(:,voxel_ind);
    data_deconv=zeros(size(bold_sig));
    event_number=zeros(1,size(bold_sig,2));
    
    disp('Retrieving HRF ...');
    if strfind(para.estimation, 'canon')
        tic
        [beta_hrf, bf, event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(bold_sig,para,temporal_mask);
        hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
        
    elseif strfind(para.estimation, 'FIR')
        tic
        para.T=1; % this needs to be = 1 for FIR
        [hrfa,event_bold] = wgr_rsHRF_FIR(bold_sig,para, temporal_mask);
        
        
    end
    
    
    
    
    nvar = size(hrfa,2); PARA = zeros(3,nvar);
    
    for voxel_id=1:nvar
        
        hrf1 = hrfa(:,voxel_id);
        
        [PARA(:,voxel_id)] = wgr_get_parameters(hrf1,para.TR/para.T);% estimate HRF parameter
        
    end
    
    toc
    disp('Done');
    
    
    
    disp('Deconvolving HRF ...');
    tic
    T = round(para.len/TR);
    if para.T>1
        hrfa_TR = resample(hrfa,1,para.T);
    else
        hrfa_TR = hrfa;
    end
    for voxel_id=1:nvar
        hrf=hrfa_TR(:,voxel_id);
        H=fft([hrf; zeros(nobs-length(hrf),1)]);
        M=fft(bold_sig(:,voxel_id));
        data_deconv(:,voxel_id) = ifft(conj(H).*M./(H.*conj(H)+2));
        event_number(voxel_id)=length(event_bold{1,voxel_id});
    end
    toc
    disp('Done');
    
    
    save(fullfile(sub_save_dir,[name,'_hrf.mat']), 'para', 'hrfa', 'event_bold', 'PARA','-v7.3');
    
    %write nifti maps of the three HRF parameters
    HRF_para_str = {'Height.nii', 'Time2peak.nii','FWHM.nii'};
    data= zeros(v.dim);
    for i=1:3
        v.fname = fullfile(sub_save_dir,[name,'_',HRF_para_str{i}]);
        data(voxel_ind)=PARA(i,:);
        spm_write_vol(v,data);
    end
    % write number of events
    v.fname = fullfile(sub_save_dir,[name,'_event_number.nii']);
    data(voxel_ind)=event_number;
    spm_write_vol(v,data);
    
    % writing back deconvolved data into nifti file
    dat3 = zeros(v1(1).dim);
    for i=1:nobs
        v1(i).fname = fullfile(sub_save_dir,[name,'_deconv',ext]);
        dat3(voxel_ind) = data_deconv(i,:);
        spm_write_vol(v1(i),dat3);
    end
end

%% example plots
figure;plot(hrfa(:,1));
figure;plot(zscore(bold_sig(:,1)));hold on;plot(zscore(data_deconv(:,1)),'r');stem(event_bold{1,1},'k');legend('BOLD','deconvolved','events')