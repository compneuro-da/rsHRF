function [beta_hrf bf event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data,xBF,temporal_mask);
%% 
% xBF.TR = 2;
% xBF.T = 8;
% xBF.T0 = fix(xBF.T/2); (reference time bin, see slice timing)
% xBF.TD_DD = 1;
% xBF.dt = xBF.TR/xBF.T;
% xBF.AR_lag = 1;
% xBF.thr = 1;
% xBF.len = 25;
% temporal_mask: generated from scrubbing.
% By: Guo-Rong Wu (gronwu@gmail.com).
% Faculty of Psychology, Southwest University.
% History:
% - 2015-04-17 - Initial version.
[N nvar] = size(data);
bf = wgr_spm_get_canonhrf(xBF);
bf2 = wgr_spm_Volterra(bf,xBF);
bf  = [bf bf2];
% xBF.name =  'Finite Impulse Response';
% xBF.length = xBF.len;
% xBF.order = 4;
%               'Gamma functions'
%               'Fourier set'
%               'Fourier set (Hanning)'
% [xBF] = spm_get_bf(xBF);
% bf = xBF.bf;
len= xBF.len;
warning off
fprintf('%d ',nvar)
% beta_hrf = sparse(nvar,xBF.TD_DD+3);
beta_hrf = cell(1,nvar);
event_bold= cell(nvar,1);
parfor i=1:nvar
%     fprintf('%d \n',i)
    [beta_hrf{1,i}  event_bold{i}] =wgr_hrf_estimation_canon(data(:,i),xBF,len,N,bf,temporal_mask);
end

beta_hrf  =cell2mat(beta_hrf);

warning on
return

function bf2 = wgr_spm_Volterra(bf,xBF)
bf2=[];
if isfield(xBF,'Volterra') 
    if xBF.Volterra==2
        bf2 = [];
        for p = 1:size(bf,2)
            for q = 1:size(bf,2)
                bf2 = [bf2 bf(:,p).*bf(:,q)];
            end
        end    
        bf2 = spm_orth(bf2);
    end
end
return

function [bf] = wgr_spm_get_canonhrf(xBF)
% Fill in basis function structure
% FORMAT [xBF] = spm_get_bf(xBF)
%
% xBF.dt      - time bin length {seconds}
% xBF.name    - description of basis functions specified
%               'hrf'
%               'hrf (with time derivative)'
%               'hrf (with time and dispersion derivatives)'
% xBF.T       - microtime resolution (for 'hrf*')

% bf      - array of basis functions
%__________________________________________________________________________
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging
% Karl Friston
% $Id: spm_get_bf.m 4473 2011-09-08 18:07:45Z guillaume $

% Edited by Guorong Wu. 2015-04-17
% Keep only double Gamma HRF.

%-Length of time bin
%--------------------------------------------------------------------------
dt = xBF.dt; %'time bin for basis functions {secs}';

%-Create basis functions
%==========================================================================
%-Microtime resolution
%----------------------------------------------------------------------
fMRI_T  =  xBF.T;
 
%-Canonical hemodynamic response function
%----------------------------------------------------------------------
[bf, p]      = spm_hrf(dt,[],fMRI_T);
p(end) = xBF.len;
[bf, p]      = spm_hrf(dt,p,fMRI_T);

%-Add time derivative
%----------------------------------------------------------------------
if xBF.TD_DD %strfind(xBF.name,'time')

    dp       = 1;
    p(6)     = p(6) + dp;
    D        = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
    bf       = [bf D(:)];
    p(6)     = p(6) - dp;

    %-Add dispersion derivative
    %------------------------------------------------------------------
    if xBF.TD_DD==2 %strfind(xBF.name,'dispersion')

        dp   = 0.01;
        p(3) = p(3) + dp;
        D    = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
        bf   = [bf D(:)];
    end
end
 
%-Orthogonalise and fill in basis function structure
%--------------------------------------------------------------------------
bf = spm_orth(bf);

return

function [beta_hrf u0]= wgr_hrf_estimation_canon(dat,xBF,len,N,bf,temporal_mask)
%% estimate HRF
thr = xBF.thr;
u0 = wgr_BOLD_event_vector(N,dat,thr,temporal_mask);
u = [    full(u0)  
    zeros(xBF.T-1,N) ];
u = reshape(u,1,[]);  %(microtime)
[beta lag] = wgr_hrf_fit(dat,len,xBF,u,N,bf);
beta_hrf = beta; beta_hrf(end+1) = lag;
return

function data = wgr_BOLD_event_vector(N,matrix,thr,temporal_mask)
%detect BOLD event.  event>thr & event<3.1
data = sparse(1,N);
if isempty(temporal_mask)
    matrix = zscore(matrix);
    for t  = 3:N-2
        if matrix(t,1) > thr && matrix(t,1) < 3.1 && all(matrix(t-2:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+2,1))% detects threshold
            data(t) = 1 ;
        end
    end
else
    datm = mean(matrix(temporal_mask));
    datstd = std(matrix(temporal_mask));
    datstd(datstd==0) = 1;%in case datstd==0;
    matrix = (matrix-datm)./datstd;
    for t  = 3:N-2
        if temporal_mask(t)
            if matrix(t,1) > thr && matrix(t,1) < 3.1 && all(matrix(t-2:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+2,1))% detects threshold
                data(t) = 1 ;
            end
        end
    end
end    

return

function  varargout = wgr_hrf_fit(dat,len,xBF,u,N,bf)
% u    - BOLD event vector (microtime).
% nlag - time lag from neural event to BOLD event .
lag = xBF.lag;
AR_lag = xBF.AR_lag;
nlag = length(lag);
erm = zeros(1,nlag);
beta = zeros(size(bf,2)+1,nlag);
% fprintf('Converged after ')
for i=1:nlag
    u_lag = [u(1,lag(i)+1:end) zeros(1,lag(i))]';
    [erm(i) beta(:,i)] = wgr_glm_estimation(dat,u_lag,bf,xBF.T,xBF.T0,AR_lag);
end
[ermin id] = min(erm); 
varargout{1} = beta(:,id);
if nargout>1
   varargout{2} = lag(id);
end
% fprintf('iterations!\n')
return

function [res_sum Beta] = wgr_glm_estimation(dat,u,bf,T,T0,AR_lag)
% u: BOLD event vector (microtime).
nscans = size(dat,1);
x = wgr_onset_design(u,bf,T,T0,nscans);
X = [x ones(nscans,1)];
[res_sum Beta] = wgr_glsco(X,dat,AR_lag);

return

function X = wgr_onset_design(u,bf,T,T0,nscans)
% u: BOLD event vector (microtime).
% bf: - basis set matrix
% T: - microtime resolution (number of time bins per scan)
% T0: - microtime onset (reference time bin, see slice timing)
ind = 1:length(u);
X     = [];
for p = 1:size(bf,2)
    x = conv(u,bf(:,p));
    x = x(ind);
    X = [X x];
end
%-Resample regressors at acquisition times
X = X( (0:(nscans - 1))*T + T0, :);

return

function [res_sum Beta] = wgr_glsco(X,Y,AR_lag)
% Linear regression when disturbance terms follow AR(p)
% -----------------------------------
% Model:
% Yt = Xt * Beta + ut , 
% ut = Phi1 * u(t-1) + ... + Phip * u(t-p) + et
% where et ~ N(0,s^2)
% -----------------------------------
% Algorithm: 
% Cochrane-Orcutt iterated regression (Feasible generalized least squares)
% -----------------------------------
% Usage:
% Y = dependent variable (n * 1 vector)
% X = regressors (n * k matrix)
% AR_lag = number of lags in AR process
% -----------------------------------
% Returns:
% Beta = estimator corresponding to the k regressors 

if nargin < 3;      AR_lag = 1;end
if nargin < 4;   max_iter = 20;end

[nobs, nvar] = size(X);
Beta = X\Y;
resid = Y - X * Beta;

max_tol = min(1e-6,max(abs(Beta))/1000);
for r = 1:max_iter    
%     fprintf('Iteration No. %d\n',r)   
    Beta_temp = Beta;

    X_AR = zeros(nobs-2*AR_lag,AR_lag);
    for m = 1:AR_lag
        X_AR(:,m) = resid(AR_lag+1-m:nobs-AR_lag-m);
    end
    
    Y_AR = resid(AR_lag+1:nobs-AR_lag);
    AR_para = X_AR\Y_AR;
    
    X_main = X(AR_lag+1:nobs,:);
    Y_main = Y(AR_lag+1:nobs);
    for m = 1:AR_lag
        X_main = X_main-AR_para(m)*X(AR_lag+1-m:nobs-m,:);
        Y_main = Y_main-AR_para(m)*Y(AR_lag+1-m:nobs-m);
    end
    
    Beta = X_main\Y_main;
    
    resid = Y(AR_lag+1:nobs) - X(AR_lag+1:nobs,:)*Beta;
    if max(abs(Beta-Beta_temp)) < max_tol
%         fprintf('%d ,',r) 
%         fprintf('Converged after %d iterations!\n',r) 
        break
    end    
    
end
res_sum = sum(resid.^2);
%if r == max_iter
%    fprintf('Maximum %d iteration reached.\n',max_iter)    
%end
return