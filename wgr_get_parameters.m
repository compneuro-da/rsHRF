function [param] = wgr_get_parameters(hdrf,dt)

% Find model parameters
%
% Height - h
% Time to peak - p (in time units of dt seconds, dt = TR/para.T)
% Width (at half peak) - w  
% Calculate Heights and Time to peak:
% modified from get_parameters2.m
param = zeros(3,1);
if any(hdrf)
    t = length(hdrf);
    n = fix(length(hdrf)*0.8);

    [h,p] = max(abs(hdrf(1:n)));
    h = hdrf(p);

    if (h >0)
        v = (hdrf >= h/2);    
    else
        v = (hdrf <= h/2);
    end;

    [a,b] = min(diff(v));
    v(b+1:end) = 0;
    w = sum(v);

    cnt = p-1;
    g =hdrf(2:end) - hdrf(1:(end-1));
    while((cnt > 0) & (abs(g(cnt)) <0.001)),
        h = hdrf(cnt);
        p = cnt;
        cnt = cnt-1;
    end;

    param(1) = h;
    param(2) = p*dt;
    param(3) = w*dt;
else
    fprintf('.')
end