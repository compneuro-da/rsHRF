clc,clear
addpath E:\rsHRF_demo_UCLA_code
%------------------------------------------------------------------------
%% White matter voxelwise HRF, basis funtions: Gamma functions
%------------------------------------------------------------------------
if 1
    % (1). HRF estimation, first unzip ./dependent_file/ds000030_csf_aCompcor5.zip into data folder ./ds000030_R1.0.5/.
    rsHRF_demo_UCLA_whitematter_HRF_batch 
    % (2). smooth results
    rsHRF_demo_UCLA_smooth_3D_wmHRF 
    % (3). Finally run 3dMVM: rsHRF_demo_UCLA_MVM_HRFpara_wm.txt
    %first copy mask file (./dependent_file/wm_mask095.nii) to /mnt/e/rsHRF_demo_UCLA/dependent_file/wm_mask095.nii
    %tcsh -x /mnt/e/rsHRF_demo_UCLA/rsHRF_demo_UCLA_MVM_HRFpara_wm.txt |& tee /mnt/e/rsHRF_demo_UCLA/out/MVM_HRFpara_wm_UCLA_diary.txt &                                                                                                              
end

%------------------------------------------------------------------------
%% Surface ROI HRF , basis funtions: Fourier set （Hanning）
%------------------------------------------------------------------------
if 1
    % (1). HRF estimation
    rsHRF_demo_UCLA_Surf_Yeo_net7_batch 
    % (2). MANOVA 
    rsHRF_demo_UCLA_Surf_ROI_manova 
end


%------------------------------------------------------------------------
%% 3D ROI FC analysis & CPM, basis funtions: Gamma functions
%------------------------------------------------------------------------
if 1
    % (1).HRF deconvolution and FC analsis
    rsHRF_demo_UCLA_DMN_ROI_FC_batch
    % (2).data preparation for prediction
    rsHRF_demo_UCLA_data_preparation4prediction
    % (3) RVR prediction and model generalization
    rsHRF_demo_UCLA_RVR_age_prediction
end
