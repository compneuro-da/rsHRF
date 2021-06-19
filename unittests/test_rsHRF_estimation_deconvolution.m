function tests = test_rsHRF_estimation_deconvolution
% Unit Tests for band filter, HRF (parameter) estimation and (iterative Wiener) deconvolution
tests = functiontests(localfunctions);

function test_rsHRF_core(testCase)
import matlab.unittest.constraints.*
fpath = fileparts(which('rsHRF.m'));
cd(fullfile(fpath,'demo_codes'))
load('HCP_100307_rfMRI_REST1_LR_Atlas_hp2000_clean_dtseries.mat')
TR = .72;
bands=[0.01 0.1]; 
bold_sig = zscore(bold_sig);
data = rsHRF_band_filter(bold_sig,TR,bands);

% Unit Tests for band filter
testCase.verifyThat(data, HasSize(size(bold_sig)));
testCase.verifyThat(data, HasLength(length(bold_sig)));

para.TR = TR;
BF = {'Canonical HRF (with time derivative)'
'Canonical HRF (with time and dispersion derivatives)'              
'Gamma functions'
'Fourier set'
'Fourier set (Hanning)'
'FIR'
'sFIR'};

para.order = 3; 
para.len = 24;

beta_size=[2 3 para.order repmat(para.order*2+1,1,2) repmat(fix(para.len/TR),1,2) ] + 2;
 
temporal_mask = []; 
para.T  = 2; 
para.T0 = 1; 
if para.T==1
    para.T0 = 1;
end

para.AR_lag = 1; 
para.thr = 1; 

min_onset_search = 4; 
max_onset_search = 8; 

for bf_id = 1:length(BF)
    para.name = BF{bf_id}; 
    if bf_id>5
        para.T  = 1; 
        para.T0 = 1; 
    end
    para.dt  = para.TR/para.T; 
    para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);

    if bf_id<=5
        [beta_hrf, bf, event_bold] = rsHRF_estimation_temporal_basis(data,para,temporal_mask); 
        hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
    else
        para.estimation = BF{bf_id}; % sFIR
        [beta_hrf, event_bold] = rsHRF_estimation_FIR(data,para,temporal_mask);
        hrfa = beta_hrf(1:end-2,:);  
    end

    
    % Unit Tests for rsHRF_estimation_temporal_basis & rsHRF_estimation_FIR
    testCase.verifyThat(beta_hrf, HasSize([beta_size(bf_id),1]));
    testCase.verifyThat(event_bold, IsOfClass('cell'));

    
    PARA = rsHRF_get_HRF_parameters(hrfa,para.dt); 
    
    % Unit Tests for rsHRF_get_HRF_parameters
    testCase.verifyThat(PARA, HasSize([3,1]));
    
    T = round(para.len/TR);
    if para.T>1
        hrfa_TR = resample(hrfa,1,para.T);
    else
        hrfa_TR = hrfa;
    end
    hrf=hrfa_TR;

    data_deconv = rsHRF_iterative_wiener_deconv(data,hrf,10,1);

    % Unit Tests for rsHRF_iterative_wiener_deconv
    testCase.verifyThat(data_deconv, HasSize(size(data)));

end