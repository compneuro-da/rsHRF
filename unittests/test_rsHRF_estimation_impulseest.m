function tests = test_rsHRF_estimation_impulseest
% Unit Tests for rsHRF_estimation_impulseest
tests = functiontests(localfunctions);

function test_rsHRF_impulseest(testCase)
import matlab.unittest.constraints.*
fpath = fileparts(which('rsHRF.m'));
cd(fullfile(fpath,'demo_codes'))
load('HCP_100307_rfMRI_REST1_LR_Atlas_hp2000_clean_dtseries.mat')
TR = .72;
bands=[0.01 0.1]; 
bold_sig = double(zscore(bold_sig));
data = rsHRF_band_filter(bold_sig,TR,bands);
para.TR = TR;
para.len = 24;
para.temporal_mask=[];
para.T  = 1; 
para.T0 = 1; 
para.thr = 1; 
min_onset_search = 4; 
max_onset_search = 8;
para.dt  = para.TR/para.T; 
para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
[beta_hrf, event_bold] = rsHRF_estimation_impulseest(data,para);
testCase.verifyThat(beta_hrf, HasSize([fix(para.len/TR)+ 1,1]));
testCase.verifyThat(event_bold, IsOfClass('cell'));