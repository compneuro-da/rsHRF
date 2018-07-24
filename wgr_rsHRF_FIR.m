function [beta_rshrf,event_bold] = wgr_rsHRF_FIR(data,para,temporal_mask)
% matlab R2015b
% temporal_mask: generated from scrubbing.
% By: Guo-Rong Wu (gronwu@gmail.com).
% Faculty of Psychology, Southwest University.
% History:
% - 2015-04-17 - Initial version.
para.temporal_mask=temporal_mask;
[N,nvar] = size(data);
if nnz(para.thr)==1
    para.thr(2)=inf;
end

% warning off
beta_rshrf = cell(1,nvar);
event_bold= cell(1,nvar);
parfor i=1:nvar
    [beta_rshrf{i}, event_bold{i}] = wgr_FIR_estimation_HRF(data(:,i),para,N);
end
beta_rshrf  =cell2mat(beta_rshrf);
% warning on
return

function [rsH,u] = wgr_FIR_estimation_HRF(data,para,N)
firmode=double(strcmp(para.estimation,'sFIR'));
if ~isfield(para,'localK')
    if para.TR<=2
        localK = 1;
    else
        localK = 2;
    end
else
    localK = para.localK;
end
u = wgr_BOLD_event_vector(N,data,para.thr,localK,para.temporal_mask);
u = find(full(u(:)));
lag = para.lag;
nlag = length(lag);
len_bin=floor(para.len/para.TR);
hrf = zeros(len_bin,nlag);
Cov_E = zeros(1,nlag);
kk=1;
for i_lag=1:nlag
    RR = u-i_lag; RR(RR<=0)=[];
    if ~isempty(RR)
        design = zeros(N,1);
        design(RR) = 1;
        [hrf(:,kk),e3] = Fit_sFIR2(data,para.TR,design,len_bin,firmode);
        Cov_E(kk) = cov(e3);
    else
        Cov_E(kk) = inf;
    end
    kk = kk+1;
end
[~, ind] = knee_pt(Cov_E); % this is a function to find the elbow point of a curve, there should be an equivalent in Python
rsH = hrf(:,ind+1);


function data = wgr_BOLD_event_vector(N,matrix,thr,k,temporal_mask)
%detect BOLD event.  event>thr & event<3.1
data = sparse(1,N);
% k=2;
if isempty(temporal_mask)
    matrix = zscore(matrix);
    for t  = 1+k:N-k
        if matrix(t,1) > thr(1) && all(matrix(t-k:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
%         if matrix(t,1) > thr(1) && matrix(t,1) < thr(2) && all(matrix(t-2:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+2,1))% detects threshold
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
            if matrix(t,1) > thr(1) && all(matrix(t-k:t-1,1)<matrix(t,1)) && all(matrix(t,1)>matrix(t+1:t+k,1))% detects threshold
                data(t) = 1 ;
            end
        end
    end
end
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
    
    b = (DX2'*DX2+sig^2*RI)\DX2'*tc;
    %     fit = DX2*b;
    e = tc - DX2*b;
    
elseif mode == 0
    
    b = pinv(DX)*tc;
    %     fit = DX*b;
    e = tc - DX*b;
    b=b(1:T);
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
