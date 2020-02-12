function [hrfa,event_bold] = rsHRF_estimation_impulseest(data,para);
% Nonparametric impulse response estimation.
% System Identification Toolbox is required.
%
% By: Guo-Rong Wu (gronwu@gmail.com).
% Faculty of Psychology, Southwest University.
% History:
% - 2015-04-17 - Initial version.
% para.thr=[1,3];
% para.lag=[2,6];
% para.TR = 2;
% para.len = 24; % second
% para.localK = 1; % local peak, f([-localK: localK]+t) <= f(t)
% if 0 % Options set for impulseest
%     options = impulseestOptions; % see impulseestOptions.m for help 
%     options.RegularizationKernel = 'HF'; %Regularizing kernel, used for regularized estimates of impulse response for all input-output channels. Regularization reduces variance of estimated model coefficients and produces a smoother response by trading variance for bias
% %                                     'TC' ? Tuned and correlated kernel, Default: 'TC'
% %                                     'none' ? No regularization is used
% %                                     'CS' ? Cubic spline kernel
% %                                     'SE' ? Squared exponential kernel
% %                                     'SS' ? Stable spline kernel
% %                                     'HF' ? High frequency stable spline kernel
% %                                     'DI' ? Diagonal kernel
% %                                     'DC' ? Diagonal and correlated kernel
%     options.PW = 5;  %Order of the input prewhitening filter. Must be one of the following:
% %                 'auto' ? Uses a filter of order 10 when RegularizationKernel is 'none'; otherwise, 0.
% %                 Nonnegative integer: Use a nonzero value of prewhitening only for unregularized estimation (RegularizationKernel is 'none').
% %                 Default: 'auto'
%     para.options = options;
% end
% [hrfa,para] = rsHRF_estimation_impulseest(randn(200,5),para);
% x = 0:para.TR:para.len;
% plot(x, hrfa)
if ~exist('impulseest.m')
   fprintf('System Identification Toolbox is required.\n') 
   return 
end
nvar = size(data,2);
if ~isfield(para,'temporal_mask') % temporal_mask: generated from scrubbing.
    para.temporal_mask=[];
end

hrfa = cell(1,nvar);
event_bold= cell(1,nvar);

parfor i=1:nvar
    [hrfa{i},event_bold{i}] = wgr_impulseest_HRF(data(:,i),para);
end
hrfa  =cell2mat(hrfa);
return

function [rsH,u] = wgr_impulseest_HRF(data,para)
nscans = size(data,1);
if ~isfield(para,'options')
    opt = [];
else
    opt = para.options;
end
if ~isfield(para,'localK')
    if para.TR>=2
        localK = 1;
    else
        localK = 2;
    end
else
    localK = para.localK;
end


u = wgr_BOLD_event_vector(nscans,data,para.thr,localK,para.temporal_mask);
u = u(:);
lag = para.lag;
nlag = length(lag);
frameRate = 1/para.TR;
NN = fix(para.len*frameRate);
bic = zeros(nlag,1);
tic
%% rs-HRF    
for i=nlag:-1:1
    u_lag = full([u(lag(i)+1:end,:); zeros(lag(i),1)]);
    srs = iddata(data, u_lag, para.TR);
%     srs = detrend(srs);
    if isempty(opt)
        srsiiu = impulseest(srs,NN);  
    else
        srsiiu = impulseest(srs,NN,opt); 
    end
    bic(i,1) = srsiiu.Report.Fit.BIC;
    HR(:,i) = srsiiu.Report.Parameters.ParVector;
end
toc
[~,id] = min(bic);
rsH = HR(:,id);
u = find(full(u)); 

function data = wgr_BOLD_event_vector(N,matrix,thr,k,temporal_mask)
%detect BOLD event.  
data = sparse(1,N);
% k = 2;
if isempty(temporal_mask)
    matrix = zscore(matrix);
    for t  = 1+k:N-k
        if matrix(t,1) > thr & all(matrix(t-k:t-1,1)<matrix(t,1)) & all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
%         if matrix(t,1) > thr & matrix(t,1) < 3.1 & all(matrix(t-k:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
            data(t) = 1 ;
        end
    end
else
    datm = mean(matrix(temporal_mask));
    datstd = std(matrix(temporal_mask));
    datstd(datstd==0) = 1;%in case datstd==0;
    matrix = (matrix-datm)./datstd;
    for t  = 1+k:N-k
        if temporal_mask(t)
            if matrix(t,1) > thr  & all(matrix(t-k:t-1,1)<matrix(t,1)) & all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
%             if matrix(t,1) > thr & matrix(t,1) < 3.1 & all(matrix(t-k:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
                data(t) = 1 ;
            end
        end
    end
end    

return
