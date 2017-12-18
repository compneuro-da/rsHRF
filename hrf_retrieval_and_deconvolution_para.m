function [sig_deconv hdrf] = hrf_retrieval_and_deconvolution(data,thr,event_lag_max,TR,T,flag)

[N nvar] = size(data);

% temporal smoothing
data(2:end-1,:) = 0.25*data(1:end-2,:) + 0.5*data(2:end-1,:) + 0.25*data(3:end,:); % you can choose whether leaving or commenting this

% raw onset.
[even_new]=wgr_trigger_onset(data,thr);


% inital HRF and its parameter.
HRF = zeros(T,nvar);
sig_deconv = single(zeros(N,nvar));
limits=round(nvar/10:nvar/10:nvar)';
limnew=zeros(length(limits)+1,1);
limnew(1)=0;limnew(2:end)=limits;
limits=limnew;clear limnew;


% adjust onset and reconstruct HRF.
warning off
if flag==1
    for ilimits=1:length(limits)-1;
        parfor i=limits(ilimits)+1:limits(ilimits+1)
            [sig_deconv(:,i) HRF(:,i) event{i} adjust_global(i)] = wgr_adjust_onset_rbeta(data(:,i),even_new{i},event_lag_max,TR,T,N);
        end
    end
    hdrf.model{1} = 'rbeta';
    hdrf.HRF{1} = HRF;
    hdrf.event{1} = event;
    hdrf.adjust_global{1} = adjust_global;
    %     hdrf.PARA{1} = PARA;
elseif flag==2
    p_m = 3; %3 - canonical + time and dispersion derivative
    for ilimits=1:length(limits)-1;
        parfor i=limits(ilimits)+1:limits(ilimits+1)
            [sig_deconv(:,i) HRF(:,i) event{i} adjust_global(i)] = wgr_adjust_onset_canonhrf(data(:,i),even_new{i},event_lag_max,TR,p_m,T,N)
        end
    end
    hdrf.model{1} = ['canonical, model type ',num2str(p_m)];
    hdrf.HRF{1} = HRF;
    hdrf.event{1} = event;
    hdrf.adjust_global{1} = adjust_global;
    %     hdrf.PARA{1} = PARA;
elseif flag==3
    mode = 1; %1 - smooth FIR
    for ilimits=1:length(limits)-1;
        parfor i=limits(ilimits)+1:limits(ilimits+1)
            [sig_deconv(:,i) HRF(:,i) event{i} adjust_global(i)] = wgr_adjust_onset_sFIR(data(:,i),even_new{i},event_lag_max,TR,mode,T,N)
        end
    end
    hdrf.model{1} = ['FIR, model type ',num2str(mode)];
    hdrf.HRF{1} = HRF;
    hdrf.event{1} = event;
    hdrf.adjust_global{1} = adjust_global;
    %     hdrf.PARA{1} = PARA;
else
    error('choose a flag value between 1 and 3');
end

warning on
clear HRF event adjust_global


% function [hrf even_new ad_global param] = wgr_adjust_onset_rbeta(dat2,even_new,event_lag_max,TR,T,N)
function [dat_deconv hrf even_new ad_global] = wgr_adjust_onset_rbeta(dat2,even_new,event_lag_max,TR,T,N)
kk=1;
hrf = zeros(T,event_lag_max+1);
Cov_E = zeros(1,event_lag_max+1);

for event_lag=0:event_lag_max
    RR = even_new-event_lag; RR(RR<=0|RR+T-1>N)=[];
    if ~isempty(RR)
        design = zeros(N,1);
        design(RR) = 1;
        matrix = repmat(RR,T,1)+repmat([0:T-1]',1,length(RR));%dat(t-past:t+future);
        hrf(:,kk) = mean(dat2(matrix),2); hrf(:,kk)=hrf(:,kk)-hrf(1,kk) ;
        s_h = conv(design,hrf(:,kk)) ;
        e3 = dat2 - s_h(1:N,1);
        Cov_E(kk) = cov(e3);
    else
        Cov_E(kk) = inf;
    end
    kk = kk+1;
end
[C ind] = min(Cov_E);
ad_global=ind-1;%begin with 0.
even_new = even_new-ad_global; even_new(even_new<=0)=[];
hrf = hrf(:,ind);
%% linear deconvolution.
H=fft([hrf; zeros(N-T,1)]);
M=fft(dat);
dat_deconv = single(ifft(conj(H).*M./(H.*conj(H)+C)));
return

% function [hrf even_new ad_global param] = wgr_adjust_onset_canonhrf(dat,even_new,event_lag_max,TR,p_m,T,N)
function [dat_deconv hrf even_new ad_global] = wgr_adjust_onset_canonhrf(dat,even_new,event_lag_max,TR,p_m,T,N)
%% global adjust.
kk=1;
hrf = zeros(T,event_lag_max+1);
for event_lag=0:event_lag_max
    RR = even_new-event_lag; RR(RR<=0)=[];
    design = zeros(N,1);
    design(RR) = 1;
    [hrf(:,kk), e3] = Fit_Canonical_HRF2(dat,TR,design,T,p_m);
    Cov_E(kk) = cov(e3);
    kk = kk+1;
end
[C ind] = min(Cov_E); ad_global=ind-1;%begin with 0.
even_new = even_new-ad_global; even_new(even_new<=0)=[];
hrf = hrf(:,ind);
%% linear deconvolution.
H=fft([hrf; zeros(N-T,1)]);
M=fft(dat);
dat_deconv = single(ifft(conj(H).*M./(H.*conj(H)+C)));
return


% function [hrf even_new ad_global param] = wgr_adjust_onset_sFIR(dat,even_new,event_lag_max,TR,mode,T,N)
function [dat_deconv hrf even_new ad_global] = wgr_adjust_onset_sFIR(dat,even_new,event_lag_max,TR,mode,T,N)
kk=1;
hrf = zeros(T,event_lag_max+1);
for event_lag=0:event_lag_max
    RR = even_new-event_lag; RR(RR<=0)=[];
    design = zeros(N,1);
    design(RR) = 1;
    [hrf(:,kk),e3] = Fit_sFIR2(dat,TR,design,T,mode);
    Cov_E(kk) = cov(e3);
    kk = kk+1;
end
[C ind] = min(Cov_E); ad_global=ind-1;%begin with 0.
even_new = even_new-ad_global; even_new(even_new<=0)=[];
hrf = hrf(:,ind);
%% linear deconvolution.
H=fft([hrf; zeros(N-T,1)]);
M=fft(dat);
dat_deconv = single(ifft(conj(H).*M./(H.*conj(H)+C)));
return



function [oneset] = wgr_trigger_onset(matrix,thr)
[N nvar] = size(matrix);
%matrix = zscore(matrix);
% Computes pseudo event.
for i = 1:nvar
    sig_vox=zscore(matrix(:,i));
    oneset_temp = [];
    for t  = 2:N-1
        if sig_vox(t) > thr && sig_vox(t-1)<sig_vox(t) && sig_vox(t)>sig_vox(t+1)% detects threshold
            oneset_temp = [oneset_temp t] ;
        end
    end
    oneset{i} = oneset_temp;
end

return


function [hrf, e] = Fit_Canonical_HRF2(tc,TR,Run,T,p)
%
% Fits GLM using canonical hrf (with option of using time and dispersion derivatives)';
%
% INPUTS:
%
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% p     - Model type
%
% Options: p=1 - only canonical HRF
%          p=2 - canonical + temporal derivative
%          p=3 - canonical + time and dispersion derivative
%
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width

len = length(Run);

X = zeros(len,p);

[h, dh, dh2] = CanonicalBasisSet(TR,T);
v = conv(Run,h);
X(:,1) = v(1:len);

if (p>1)
    v = conv(Run,dh);
    X(:,2) = v(1:len);
end

if (p>2)
    v = conv(Run,dh2);
    X(:,3) = v(1:len);
end

X = [(zeros(len,1)+1) X];
b = pinv(X)*tc;
e = tc-X*b;
% fit = X*b;

b = b(2:end);

if (p == 2)
    bc = sign(b(1))*sqrt(b(1)^2 + b(2)^2);
    H = [h dh];
elseif (p==1)
    bc = b(1);
    H = h;
elseif (p>2)
    bc = sign(b(1))*sqrt(b(1)^2 + b(2)^2 + b(3)^2);
    H = [h dh dh2];
end

hrf = H*b;


return
% END MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, dh, dh2] = CanonicalBasisSet(TR,len)

% len = round(30/TR);
xBF.dt = TR;
xBF.length= len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);

v1 = xBF.bf(1:len,1);
v2 = xBF.bf(1:len,2);
v3 = xBF.bf(1:len,3);

h = v1;
dh =  v2 - (v2'*v1/norm(v1)^2).*v1;
dh2 =  v3 - (v3'*v1/norm(v1)^2).*v1 - (v3'*dh/norm(dh)^2).*dh;

h = h./max(h);
dh = dh./max(dh);
dh2 = dh2./max(dh2);

return

function [hrf, e] = Fit_sFIR2(tc,TR,Runs,T,mode)
% function [hrf, fit, e, param] = Fit_sFIR(tc,TR,Runs,T,mode)
%
% Fits FIR and smooth FIR model
%
% INPUTS:
%
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% mode  - FIR or smooth FIR
%   options:
%       0 - standard FIR
%       1 - smooth FIR
%
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 10/02/09

[DX] = tor_make_deconv_mtx3(Runs,T,1);
DX2 = DX(:,1:T);
num = T;

if mode == 1
    
    C=(1:num)'*(ones(1,num));
    h = sqrt(1/(7/TR));                       % 7 seconds smoothing - ref. Goutte
    
    v = 0.1;
    sig = 1;
    
    R = v*exp(-h/2*(C-C').^2);
    RI = inv(R);
    
    b = inv(DX2'*DX2+sig^2*RI)*DX2'*tc;
    %     fit = DX2*b;
    e = tc - DX2*b;
    
elseif mode == 0
    
    b = pinv(DX)*tc;
    %     fit = DX*b;
    e = tc - DX*b;
    
end

hrf = b;

return


function [DX,sf] = tor_make_deconv_mtx3(sf,tp,eres,varargin)
% function [DX,sf] = tor_make_deconv_mtx(sf,tp,eres,[opt] TRs before stim onset,[num. sessions],[docenter],[scanspersess])
%   sf: cell array of stick functions, one per condition
%       all sf cells should be of the same length
%	Or matrix of stick functions, 1 column per condition
%
%
%   tp: number of timepoints to estimate in hrf deconvolution matrix
%   eres: timebins in sf array for each TR
%
%   DX: deconvolution matrix
%       estimates O.tp time points for each condition
%       Time resolution is in TRs
%
%   sf: stick function resampled at TR
%
%   Optional:
%   1 - TRs before: 0 or number of time-points to shift LEFT
%   2 - number of sessions; if > 1, adds session-specific intercepts
%   3 - docenter, 1/0 for do/do not center columns, default 0
%   4 - scanspersess: how many scans per session?  prevents regressors from
%       running over into the next session (recursive).
%
%
%   No parametric modulation of sf's allowed.
%
% Tor Wager, 10/20/01   modified 9/20/02 for variable tp's for diff evt types
% modified 4/22/04  to center columns, 2/9/05 for multi-session boundary
% respect

docenter = 0;
if nargin > 5, docenter = varargin{3};,end

if ~iscell(sf)
    %sf = mat2cell(sf,size(sf,1),ones(size(sf,2)));
    for i = 1:size(sf,2), sf2{i} = sf(:,i);, end
    sf = sf2;
end

if length(tp) == 1, tp = repmat(tp,1,length(sf));, end
if length(tp) ~= length(sf), error('timepoints vectors (tp) and stick function (sf) lengths do not match!'), end

tbefore = 0;
nsess = size(sf,1);

if nargin > 4, nsess = varargin{2};, end
if nargin > 3, tbefore = varargin{1};, end

shiftElements = eres;
% each time point is timeRes TRs.


% multiple sessions
% prevent regressors from running across session lines
if nargin > 6, numframes = varargin{4};,
    
    st = cumsum([1 numframes]);
    en = st(2:end) - 1;         % ending values
    st = st(1:end-1);           % starting values
    
    for sess = 1:length(numframes)
        
        % get ons for this session only
        for i = 1:length(sf)
            sfsess{i} = sf{i}(st(sess):en(sess));
        end
        
        % get DX for this session only, no centering
        [DXs{sess,1}] = tor_make_deconv_mtx3(sfsess,tp,eres,varargin{1},varargin{2},0);
    end
    
    % concatenate across sessions
    DX = cat(1,DXs{:});
    
    % intercepts
    %warning('intercept term removed');
    DX = DX(:,1:end-1);                     % remove overall intercept
    DX = [DX intercept_model(repmat(numframes, 1, length(numframes)))];   % get session-specific
    
else
    % run the single-session model
    
    
    % -------------------------------------------------------------------
    % * downsample sf to number of TRs
    % -------------------------------------------------------------------
    numtrs = round(length(sf{1}) ./ eres);
    myzeros = zeros(numtrs,1);
    origsf = sf;
    
    for i = 1:length(sf)
        Snumtrs = length(sf{i}) ./ eres;
        if Snumtrs ~= round(Snumtrs), warning(['sf{ ' num2str(i) '}: length not evenly divisible by eres.']),end
        if numtrs ~= Snumtrs, warning(['sf{ ' num2str(i) '}: different length than sf{1}.']),end
        
        inums = find(sf{i} > 0);
        inums = inums ./ eres;      % convert to TRs
        inums = ceil(inums);        % nearest TR
        % i'm getting weird effects with round
        % and stimuli onsets in between TRs -
        % if the 1 in the matrix occurs before the
        % onset of the actual event, the est. HRF is inverted -
        % thus i'm now using ceiling.
        inums(inums == 0) = 1;      % never use 0th element
        sf{i} = myzeros;
        sf{i}(inums) = 1;           % always use 1 for sf
    end
    
    % plot to check sampling of delta function
    %figure;
    %for i = 1:length(sf)
    %    subplot(length(sf),1,i);hold on
    %    plot(1:1/length(origsf{i}):2-1/length(origsf{i}),origsf{i})
    %    plot(1:1/length(sf{i}):2-1/length(sf{i}),sf{i},'r')
    %end
    
    % -------------------------------------------------------------------
    % * make deconvolution matrix DX
    % -------------------------------------------------------------------
    
    index = 1;
    for i = 1:size(sf,2)
        
        if tbefore ~= 0
            for j = tbefore:-1:1
                mysf = [sf{i}(j+1:end); zeros(j,1)];
                DX(:,index) = mysf;
                index = index + 1;
            end
        end
        
        DX(:,index) = sf{i};
        index = index + 1;
        inums = find(sf{i} == 1);
        
        for j = 2:tp(i)
            inums = inums + 1;      % + 1 because we've downsampled already.  + shiftElements;
            reg = myzeros;
            reg(inums) = 1;
            reg = reg(1:numtrs);
            while length(reg) < size(DX,1), reg = [reg;0];,end % add 0's if too short
            try
                DX(:,index) = reg;
            catch
                whos DX
                whos reg
                error('Different column lengths!')
            end
            index  = index + 1;
        end
        
    end
    
    % -------------------------------------------------------------------
    % * add intercept
    % -------------------------------------------------------------------
    if nsess < 2
        DX(:,end+1) = 1;
    else
        index = 1;
        scanlen = size(DX,1) ./ nsess;
        if round(scanlen) ~= scanlen, warning('Model length is not an even multiple of scan length.'),end
        
        for startimg = 1:scanlen:size(DX,1)
            X(startimg:startimg+scanlen-1,index) = 1;
            index = index + 1;
        end
        
        DX = [DX X];
    end
    
    
end     % multi-session vs. single

if docenter
    % center columns (not intercepts)
    wh = 1:size(DX,2)-nsess;
    DX(:,wh) = DX(:,wh) - repmat(mean(DX(:,wh)),size(DX,1),1);
end


return

