function para = wgr_rsHRF_global_para()
%% miminum fMRI volumes (3D)
para.volume_threshold = 30;

%% Number of basis functions (k),Only setting for ''Gamma Functions'' (k), ''Fourier Set'' (2k+1), ''Fourier Set (Hanning)'' (2k+1)';
para.num_basis = 4;

%% save pvalue for seed based GC map (NIfTI/GIfTI)
para.flag_pval_pwgc = 1;

para.regmode = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR')

%% save response height PSC (percent signal change) maps (NIfTI/GIfTI)
para.flag_save_psc = 1;

%% Combine all input signals for connectivity analysis
para.combine_ROI = 1; % 1 combine, 0: no

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
