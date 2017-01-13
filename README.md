Matlab code for resting state HRF estimation
========
![BOLD HRF](https://github.com/guorongwu/rsHRF/raw/master/docs/BOLD_HRF.png)

**PLEASE visit <https://guorongwu.github.io/HRF/>  for detail information for resting-state HRF deconvolution.**

Quickstart 
-------------
(canon2dd: canonical HRF with its delay and dispersion derivatives) 

**BOLD fMRI** parameters setting
```
temporal_mask = []; % without mask, it means temporal_mask = ones(nobs,1); i.e. all time points included. nobs: number of observation = size(data,1). if want to exclude the first 1~5 time points, let temporal_mask(1:5)=0;
```
```
data: nobs x nvar (nvar: number of variables; e.g. 200x90, 200x 50000, ....)
```
```
TR = 2;

para.TR = TR;

para.T  = 5; % temporal grid: TR/5. magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid

para.T0 = 3; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)

para.dt  = para.TR/para.T; % fine scale time resolution.

para.TD_DD = 2; % time and dispersion derivative

para.AR_lag = 1; % AR(1) noise autocorrelation.

para.thr = 1; % (mean+) para.thr*standard deviation threshold to detect event.

para.len = 24; % length of HRF, here 24 seconds

para.lag  = fix(3/para.dt):fix(9/para.dt); % 3 to 9 seconds
```
HRF estimation

```
[beta_hrf bf event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data,para,temporal_mask);
```
```
hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
```
HRF parameters estimation (PARA)

```
hrf1 = hrfa(:,1); 

plot(hrf1) % HRF shape visualisation

% do a for loop for other variable: 

nvar = size(hrfa,2); PARA = zeros(3,nvar);

for i=1:nvar; 

	hrf1 = hrfa(:,i); 
	
	[PARA(:,i)] = wgr_get_parameters(hrf1,para.TR/para.T);% estimate HRF parameter 
	
end
```
```
PARA(1,:): response height (response magnitude of neuronal activity)
PARA(2,:): Time to peak (latency of neuronal activity)
PARA(3,:): Width / FWHM (duration of neuronal activity)
Response height (percent signal change) = PARA(1,:)./beta_hrf(end-1,:)*100; 

```



**Citation**
--------

_Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding, Huafu Chen, Daniele Marinazzo*. "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data." Medical Image Analysis, 2013, 17:365-374. [PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/2013_MIA.pdf)_

_Guo-Rong Wu, Daniele Marinazzo. "Sensitivity of the resting state hemodynamic response function estimation to autonomic nervous system fluctuations." Philosophical Transactions of the Royal Society A, 2016, 374: 20150190.[PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/2016_PTA.pdf)_

_Guo-Rong Wu, Daniele Marinazzo. "Retrieving the Hemodynamic Response Function in resting state fMRI: methodology and applications." PeerJ PrePrints, 2015.[PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/poster_OHBM2016_HRF.pdf)_
