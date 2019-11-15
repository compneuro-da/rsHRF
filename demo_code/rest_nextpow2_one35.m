function Result = rest_nextpow2_one35(n)
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%------------------------------------------------------------------------------------------------------------------------------
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song 
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

    if length(n)>1
        n = cast(length(n),class(n));
    end
    if n<16
        Result =2^nextpow2(n);
        return;
    end 
    
    limit =nextpow2(n);             %n=134, limit=8
    tbl=[2^(limit-1):2^limit];      %tbl =128, 129, ... , 256
    tbl =tbl(find(tbl>=n));          %tbl =134, 135, ... , 256
    for x=1:length(tbl)
        Result =tbl(x);
        [f,p]=log2(Result);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
        if mod(Result,3*5)==0        
            y= Result /(3*5);
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,3)==0        
            y= Result /3;
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,5)==0        
            y= Result /5;
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
    end
    Result =NaN;    % Should not reach, except when n=1

% csfft_nextup35 in AFNI list 1~1024, 20070516, dawnsong
% 2
% 4
% 6
% 8
% 10
% 12
% 16
% 20
% 24
% 30
% 32
% 40
% 48
% 60
% 64
% 80
% 96
% 120
% 128
% 160
% 192
% 240
% 256
% 320
% 384
% 480
% 512
% 640
% 768
% 960
% 1024