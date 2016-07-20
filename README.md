Matlab code for resting state HRF estimation
========

Demo
----

TR = 2;

para.TR = TR;

para.T = 11;%

para.T0 = 6; % slice ref time, I always reference to middle slice time (preprocessing: slice timing)

para.dt     = para.TR/para.T; % fine scale time resolution.

para.TD_DD = 2; % time and dispersion derivative

para.AR_lag = 2; % AR(2) noise autocorrelation.

para.thr = 1; % SD threshold to detect event.

para.len = 24; % length of HRF, here 24 seconds

para.lag  = fix(3/para.dt):fix(9/para.dt); % 3 to 9 seconds

temporal_mask = []; % without mask, it's equal to temporal_mask = ones(nobs,1); nobs: number of observation = size(data,1).
% if want to exclude first 1~5 time points let temporal_mask(1:5)=0;

‘’’
[beta_hrf bf event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data,para,temporal_mask);
‘’’

hrf = bf*beta_hrf(1:size(bf,2),:); %HRF

call the code to calculate HRF parameters: PARA
response height (percent signal change) = PARA(1)./beta_hrf(end-1,:)*100; 
this is more meaningful than PARA(1), the map on brain is not like the one we see in MIA paper.
but you can also use PARA(1)---it's a absolute value. 



**Citation**
--------
. Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding, Huafu Chen, Daniele Marinazzo*. "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data." Medical Image Analysis, 2013, 17:365-374.

. Guo-Rong Wu, Daniele Marinazzo. "Sensitivity of the resting state hemodynamic response function estimation to autonomic nervous system fluctuations." Philosophical Transactions of the Royal Society A, 2016, 374: 20150190.

. Guo-Rong Wu, Daniele Marinazzo. "Retrieving the Hemodynamic Response Function in resting state fMRI: methodology and applications." PeerJ PrePrints, 2015.