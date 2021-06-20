clc,clear
%CIfTI toolbox
addpath E:\HCP_Retest\cifti-matlab-master
%SimTB toolbox
addpath E:\simtb_v18\sim
%data from HCP
dat=cifti_read('E:\HCP_Retest\103818\MNINonLinear\Results\rfMRI_REST1_LR\rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
data = dat.cdata';
TR = 0.72;
dataf = rsHRF_band_filter(data,TR,[0.01 0.1]);
lag = fix(6/TR);

xBF.TR = TR;
xBF.T = 1;
xBF.T0 = 1;
xBF.dt = xBF.TR/xBF.T;
xBF.AR_lag = 0;
xBF.thr = 1;
xBF.len = 25;
xBF.localK = 2;
xBF.order = 3;
xBF.name = 'Gamma functions';
xBF.lag = fix(4/xBF.dt):fix(9/xBF.dt);

data_length = 40:40:1200;   % data length
nv = 1:length(data_length); % number of events
snr =  0;
nkk = 2000;
for hh=1:length(nv)
    k=1; beta_hrf={};
    for i=1:size(dataf,2)
        evBOLD{i}=rsHRF_find_event_vector(dataf(1:data_length(hh),i),1,2,[]);
        tmp = find(evBOLD{i})-lag; tmp(tmp<=0)=[];
        tmp2 = zeros(data_length(hh),1); tmp2(tmp)=1;
        evNeu{i} = tmp2;
        numev(i,1) = length(tmp);
        if numev(i,1)==nv(hh) % 
            [TC, MDESC, P] = simtb_TCsource(tmp2, TR, 1); % all use same HRF
            TCn = awgn(TC,snr,'measured');
            TCnf = rsHRF_band_filter(TCn,TR,[0.01 0.1]);
            evB = rsHRF_find_event_vector(TCnf,1,2,[]);
            if nnz(evB)==nv(hh)
                evNeu4{k,1}= tmp2;
                [beta_hrf{k,1}, bf] = rsHRF_estimation_temporal_basis(TC,xBF,[],0);
                k=k+1;
                if k>nkk%2000
                    break;
                end
            end
        end
    end
    hrf = simtb_spm_hrf(TR, P);
    hrfa = cell2mat( arrayfun(@(x) bf*beta_hrf{x}(1:end-2),1:length(beta_hrf),'UniformOutput', false) );

    nobs = size(hrfa,1);
    xx = (1:nobs)*xBF.dt;
    mf0 = mean(hrfa,2);
    sf0 = std(hrfa,0,2);
    
    nk(hh) = k;
    hrf_all(:,hh) = hrf;
    hrfm_all(:,hh) = mf0;
    hrfs_all(:,hh) = sf0;
    mf(hh,:)  = [mean(mf0) mean(hrf)];
    cv(hh,:)  = [mean(sf0)/mean(mf0) mean(sf0)/mean(hrf)];

end
save(['sim_hrf_HCP_SNR',num2str(snr),'_',num2str(nkk),'trial','.mat'],'hrf*_all','mf','cv');

figure('color','w');plot(cv);
xlabel({['# events (Add white Gaussian noise to signal,SNR=',num2str(snr),')'],'TR=0.72s, Data Length (40:40:1200~ #events 1:30)'})
ylabel('Coefficient of Variation (sigma/mu)')
for i=1:2
    [~,mmm]=rsHRF_knee_pt(cv(:,i));
    hold on;plot(mmm,cv(mmm,i),'o')
end
legend({'mu = mean(estimated HRF)','mu = mean(ground truth HRF)'})
