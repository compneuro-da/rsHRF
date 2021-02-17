function output = prt_machine_rvr(d,args)
% Relevance vector regression (training and testing)
% FORMAT output = prt_machine_svm_bin(d,args)
% Inputs:
%   d         - structure with data information, with mandatory fields:
%     .train      - training data (cell array of matrices of row vectors,
%                   each [Ntr x D]). each matrix contains one representation
%                   of the data. This is useful for approaches such as
%                   multiple kernel learning.
%     .test       - testing data  (cell array of matrices row vectors, each
%                   [Nte x D])
%     .tr_targets - training labels (for classification) or values (for
%                   regression) (column vector, [Ntr x 1])
%     .use_kernel - flag, is data in form of kernel matrices (true) of in
%                form of features (false)
%    args     - libSVM arguments
% Output:
%    output  - output of machine (struct).
%     * Mandatory fields:
%      .predictions - predictions of classification or regression [Nte x D]
%     * Optional fields:
%      .func_val - value of the decision function
%      .type     - which type of machine this is (here, 'classifier')
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Carlton Chu
% $Id: prt_machine_rvr.m 741 2013-07-22 14:07:36Z mjrosa $


output.type=d.pred_type;

SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a string (empty or otherwise)

    K=d.train{1};
    [n m]=size(K);
    if n~=m
        error('prt_machine_krr:kernelSize',['Error: krr'...
            ' training kernel should be square ' ...
            ' SOLUTION: do the right thing']);
    end
    
    % Run RVR
    %--------------------------------------------------------------------------
    t=d.tr_targets;
    w=prt_rvr(K,t);
    
    output.predictions=d.test{1}*w(1:end-1)+w(end);
    output.func_val=output.predictions;
    output.alpha=w(1:end-1);
  
end

end
