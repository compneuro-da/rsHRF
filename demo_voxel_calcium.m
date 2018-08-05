%% demo code for voxel-wise HRF deconvolution
%% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
%% Guo-Rong Wu, gronwu@gmail.com, UESTC, UGent, 2013.9.12
%% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D..
%% A blind deconvolution approach to recover effective connectivity brain networks
%% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .
clc,clear;warning off all
load calcium
bold_sig=BOLD_calcium';
%%===========PARAMETERS========================
% choose the set of basis functions THIS MUST BE AN INPUT

%para.estimation='canon2dd'; % this one for canonical HRF plus two derivatives
para.estimation='sFIR'; % this one for smoothed FIR
%para.estimation='FIR'; % this one for unsmoothed FIR

temporal_mask = []; % without mask, it means temporal_mask = ones(nobs,1); i.e. all time points included. nobs: number of observation = size(data,1). if want to exclude the first 1~5 time points, let temporal_mask(1:5)=0;

TR = 1.5; % THIS WILL BE READ FROM THE BIDS DATA

para.TR = TR;

para.passband=[0.01 0.08]; %bandpass filter lower and upper bound

%%% the following parameter (upsample grid) can be > 1 only for Canonical.
%%% Set = 1 for FIR
para.T  = 1; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid

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
nobs=size(bold_sig,1);
bold_sig=zscore(bold_sig);
bold_sig = rest_IdealFilter(bold_sig, para.TR, para.passband);
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
hrf=hrfa_TR;
H=fft([hrf; zeros(nobs-length(hrf),1)]);
M=fft(bold_sig(:,1));
data_deconv = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));
event_number=length(event_bold{1,1});
toc
disp('Done');



%% example plots
event_plot=sparse(1,nobs);
event_plot(event_bold{1,1})=1;
figure;plot((1:length(hrfa(:,1)))*TR,hrfa(:,1));xlabel('time (s)')
figure;plot(BOLD_time,zscore(bold_sig(:,1)));
hold on;plot(calcium_time,zscore(calcium_raw)-10,'g')
%plot(calcium_time,trigger_times-3,'m')
plot(BOLD_time,zscore(data_deconv(:,1)),'r');
stem(BOLD_time,event_plot,'k');
xlabel('time (s)')
legend('BOLD','calcium','deconvolved','events')