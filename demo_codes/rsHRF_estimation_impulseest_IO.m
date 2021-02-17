function [hrf] = rsHRF_estimation_impulseest_IO(input, output, para)
%% detrend
if ~isempty(para.flag_detrend)
    input = spm_detrend(input,para.flag_detrend);
    output = spm_detrend(output,para.flag_detrend);
end
%% filtering
if ~isempty(para.Band)
    output = bandpass_filt(output, 1/para.TR, para.Band, 10, 'fir', 'twopass', 1);
%     output = rsHRF_band_filter(output,para.TR,para.Band);
end
 
if 0
    figure
    plot(zscore([input output]))
end

%% neural binay stimulus onset
u0 = rsHRF_find_event_vector(input,para.thr,para.localK,para.temporal_mask);
u0 = full(u0(:));
srs0 = iddata(output, u0, para.TR);
srs0.InputName  = 'LFP-bin';
srs0.OutputName = 'OIS';
NN=fix(para.len*1/para.TR);
srsiiu_lfp_bin = impulseest(srs0,NN,para.options);
 
%% LFP HRF
srs = iddata(output, input, para.TR);
srs.InputName  = 'LFP';
srs.OutputName = 'OIS';
srsii_lfp = impulseest(srs,NN,para.options);

H_LFP_bin = srsiiu_lfp_bin.Report.Parameters.ParVector;
H_LFP = srsii_lfp.Report.Parameters.ParVector;

[H_Blind] = rsHRF_estimation_impulseest(output,para);
hrf.H_LFP_bin=H_LFP_bin;
hrf.H_LFP=H_LFP;
hrf.H_Blind=H_Blind;




function filt = bandpass_filt(data, Fs, Fbp, N, type, dir, meanback)
% 
% BANDPASS_FILT applies a band-pass filter to the data in a specified
% frequency band.
% 
% 
%USAGE
%-----
% filt = bandpass_filt(data, Fsample, Fbp, N, type, dir, meanback)
% 
% 
%INPUT
%-----
% - DATA    : data matrix (Nsamples x Nchannels)
% - FSAMPLE : sampling frequency in Hz
% - FBP     : frequency band, specified as [Fhp Flp]
%             Fhp = 0   => low-pass filter
%             Flp = inf => high-pass filter
% - N       : filter order (positive integer). If N=[], the default value
%             will be used: 4 ('but') or dependent upon frequency band and
%             data length ('fir', 'firls')
% - TYPE    : filter type. If TYPE='', the default value will be used.
%             'but'   Butterworth IIR filter (default)
%             'fir'   FIR filter using Matlab fir1 function
%             'firls' FIR filter using Matlab firls function (requires
%              Matlab Signal Processing Toolbox)
% - DIR     : filter direction. If DIR='', the default value will be used.
%             'onepass'         forward filter only
%             'onepass-reverse' reverse filter only, i.e. backward in time
%             'twopass'         zero-phase forward and reverse filter (default)
%             'twopass-reverse' zero-phase reverse and forward filter
%             'twopass-average' average of the twopass and the twopass-reverse
% - MEANBACK: 1 (add the mean back after filtering [default]) or 0 (do not add)
%
% Note that a one- or two-pass filter has consequences for the strength of
% the filter, i.e. a two-pass filter with the same filter order will
% attenuate the signal twice as strong.
% 
% 
%EXAMPLE
%-------
% T  = .5; % [s] => sampling freq. = 1/T => Nyquist freq. = 1/(2T)
% t  = (0:T:2*pi*10)';
% f1 = (1/2/T)*.030;
% f2 = (1/2/T)*.050;
% f3 = (1/2/T)*.075;
% x1 = 10*sin(2*pi*f1*t);
% x2 = 4*sin(2*pi*f2*t);
% x3 = 15*cos(2*pi*f3*t);
% x = x1 + x2 + x3;
% 
% type = {'but' 'fir' 'firls'};
% dir = {'onepass' 'onepass-reverse' 'twopass' 'twopass-reverse' 'twopass-average'};
% addmean = 0; N=10; tt=1; dd=3;
% 
% xLP = bandpass_filt(x, 1/T, [0 (f2+f3)/2], N, type{tt}, dir{dd}, addmean);
% figure, plot(t,x,t,xLP,t,x1+x2)
% title(sprintf('Low-pass filter (%s - %s)', type{tt}, dir{dd}))
% legend('Original', 'Filtered', 'f1 & f2')
% xHP = bandpass_filt(x, 1/T, [(f1+f2)/2 inf], N, type{tt}, dir{dd}, addmean);
% figure, plot(t,x,t,xHP,t,x2+x3)
% title(sprintf('High-pass filter (%s - %s)', type{tt}, dir{dd}))
% legend('Original', 'Filtered', 'f2 & f3')
% xBP = bandpass_filt(x, 1/T, [(f1+f2)/2 (f2+f3)/2], 6, type{tt}, dir{dd}, addmean);
% figure, plot(t,x,t,xBP,t,x2)
% title(sprintf('Band-pass filter (%s - %s)', type{tt}, dir{dd}))
% legend('Original', 'Filtered', 'f2')

% Copyright (c) 2003-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_preproc_lowpassfilter.m 7123 2012-12-06 21:21:38Z roboos $
% $Id: ft_preproc_bandpassfilter.m 7123 2012-12-06 21:21:38Z roboos $
% $Id: ft_preproc_highpassfilter.m 7123 2012-12-06 21:21:38Z roboos $
% $Id: filter_with_correction.m 7123 2012-12-06 21:21:38Z roboos $


% Default values
%==========================================================================
if nargin<4 || isempty(N)
    N = []; % set the default filter order later
end
if nargin<5 || isempty(type)
    type = 'but'; % set the default filter type
end
if nargin<6 || isempty(dir)
    dir = 'twopass'; % set the default filter direction
end
if nargin<7 || isempty(meanback)
    meanback = 1; % set the default option for the mean
end


% Initialize
%==========================================================================
Fbp = sort(Fbp);
if Fbp(1)==0 && isinf(Fbp(2)) % full spectrum => nothing should be done
    filt = data;
    return
end
if Fbp(1)==0 % low-pass filter
    Fbp   = Fbp(2);
    ftype = -1; % filter type
elseif isinf(Fbp(2)) % high-pass filter
    Fbp   = Fbp(1);
    ftype = 1; % filter type
else % band-pass filter
    ftype = 0; % filter type
end
Fn = Fs/2; % Nyquist frequency
Nsamples = size(data, 1);
data = data.'; % Nchannels x Nsamples


% Compute filter coefficients
%==========================================================================
switch type
    case 'but'
        if isempty(N)
            if ftype==-1 || ftype==1 % low- and high-pass filters
                N = 6;
            elseif ftype==0 % band-pass filter
                N = 4;
            end
        end
        if ftype==-1 % low-pass filter
            [B, A] = butter(N, Fbp/Fn);
        elseif ftype==1 % high-pass filter
            [B, A] = butter(N, Fbp/Fn, 'high');
        elseif ftype==0 % band-pass filter
            [B, A] = butter(N, [Fbp(1)/Fn Fbp(2)/Fn]);
        end
    case 'fir'
        if isempty(N)
            if ftype==-1 || ftype==1 % low- and high-pass filters
                N = 3*fix(Fs / Fbp);
                if ftype==1 && rem(N, 2)==1
                    N = N + 1;
                end
            elseif ftype==0 % band-pass filter
                N = 3*fix(Fs / Fbp(1)); % old: 25
            end
        end
        if N > floor( (Nsamples - 1) / 3)
            if ftype==1
                N = floor(Nsamples/3) - 2;
                if rem(N, 2)==1
                    N = N + 1; % "fir1" always uses an even filter order for the highpass configuration
                end
            else
                N = floor(Nsamples/3) - 1;
            end
        end
        if ftype==-1 % low-pass filter
            [B, A] = fir1(N, Fbp/Fn);
        elseif ftype==1 % high-pass filter
            [B, A] = fir1(N, Fbp/Fn, 'high');
        elseif ftype==0 % band-pass filter
            [B, A] = fir1(N, [Fbp(1)/Fn Fbp(2)/Fn]);
        end
    case 'firls' % from NUTMEG's implementation
        if isempty(N)
            if ftype==-1 || ftype==1 % low- and high-pass filters
                N = 3*fix(Fs / Fbp);
                if ftype==1 && rem(N, 2)==1
                    N = N + 1;
                end
            elseif ftype==0 % band-pass filter
                N = 3*fix(Fs / Fbp(1));
            end
        end
        if N > floor( (Nsamples - 1) / 3)
            if ftype==1
                N = floor(Nsamples/3) - 2;
                if rem(N, 2)==1
                    N = N + 1;
                end
            else
                N = floor(Nsamples/3) - 1;
            end
        end
        f = 0:0.001:1;
        if rem(length(f), 2)~=0
            f(end) = [];
        end
        z = zeros(1, length(f));
        if ftype==-1 % low-pass filter
            [val, pos1] = min(abs(Fs*f/2 - 0));
            if isfinite(Fbp)
                [val, pos2] = min(abs(Fs*f/2 - Fbp));
            else
                pos2 = length(f);
            end
        elseif ftype==1 % high-pass filter
            [val,pos1] = min(abs(Fs*f/2 - Fbp));
            pos2 = length(f);
        elseif ftype==0 % band-pass filter
            if isfinite(Fbp(1))
                [val, pos1] = min(abs(Fs*f/2 - Fbp(1)));
            else
                [val, pos2] = min(abs(Fs*f/2 - Fbp(2)));
                pos1 = pos2;
            end
            if isfinite(Fbp(2))
                [val, pos2] = min(abs(Fs*f/2 - Fbp(2)));
            else
                pos2 = length(f);
            end
        end        
        z(pos1:pos2) = 1;
        A = 1;
        B = firls(N, f, z); % requires Matlab signal processing toolbox
    otherwise
        error('Unsupported filter type "%s"', type);
end


% Apply A to the data and correct edge-artifacts for one-pass filtering
%==========================================================================
poles = roots(A);
if any(abs(poles) >= 1)
    %poles
    fprintf(['Calculated filter coefficients have poles on or outside\n' ...
        'the unit circle and will not be stable. Try a higher cutoff\n' ...
        'frequency or a different type/order of filter.\n'])
    filt = [];
    return
end

dcGain = sum(B)/sum(A);
mu     = sum(data, 2)/Nsamples;    % data mean
data   = bsxfun(@minus, data, mu); % zero-mean

switch dir
    case 'onepass'
        offset = data(:,1);
        data   = data - repmat(offset, 1, Nsamples);
        filt   = filter(B, A, data.').' + repmat(dcGain*offset, 1, Nsamples);
        % old: filt = filter(B, A, data.');
    case 'onepass-reverse'
        offset = data(:, end);
        data   = fliplr(data) - repmat(offset, 1, Nsamples);
        filt   = filter(B, A, data.').';
        filt   = fliplr(filt) + repmat(dcGain*offset, 1, Nsamples);
        % old: data=fliplr(data); filt=filter(B, A, data.'); filt fliplr(filt);
    case 'twopass'
        % filtfilt does the correction for us
        filt = filtfilt(B, A, data.').';
    case 'twopass-reverse'
        % filtfilt does the correction for us
        filt = fliplr(filtfilt(B, A, fliplr(data).').');
    case 'twopass-average'
        % take the average from the twopass and the twopass-reverse
        filt1 = filtfilt(B, A, data.').';
        filt2 = fliplr(filtfilt(B, A, fliplr(data).').');
        filt  = (filt1 + filt2)/2;
    otherwise
        error('Unsupported filter direction "%s"', dir);
end


if meanback
    filt = bsxfun(@plus, filt, mu); % add mean back to the filtered data
end

filt = filt.'; % restore data shape