function para = wgr_rsHRF_global_para()
% denoise
if 0
    % global regression
    para.globa_reg = 1;
    para.aCompcor = 0;
else
    % acompcor
    para.globa_reg = 0;
    para.aCompcor = 1;
    para.aCompcor_numcomps = 5;
end

%%  the significance level for determination of outliers
para.pvalue_rm = 0.05;
para.rmoutlier_deconv = 1; %replace outlier in deconvolved signals (only work for 4D data)
para.Inpainted = 0; %Inpainted NAN, using the plate method
% para.Inpainted = 1; %Inpainted NAN with the springs method. It is often considerably faster, but not quite as smooth.

%% local spontaneou event 
para.localK = 2; % local peak, f([-localK: localK]+t) <= f(t)

%% whether delete temporary files (brain mask, atlas in fMRI spaces)
para.delete_files = 1;

%% default mask: for human brain (in normalizedd MNI space)
isMNI = 0;
if isMNI
    para.mask_nii = fullfile(spm('Dir'),'tpm','mask_ICV.nii');
    if ~exist(para.mask_nii,'file')
        para.mask_nii ='';
    end
else
    para.mask_nii ='';
end
%% default mask: animal?
isanimal = 0; % animal ?
if isanimal
    para.mask_nii ='';
end