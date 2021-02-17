clc,clear
main = 'E:\data\ds000030_R1.0.5\derivatives\fmriprep';
mainoutdir = 'E:\data\UCLA_ROI9_FC';
MNI=[0,-52,7;
    -1,54,27;
    -46,-66,30;
    49,-63,33;
    -61,-24,-9;
    58,-24,-9;
    0,-12,9;
    -25,-81,-33;
    25,-81,-33];
ROI = [MNI ones(size(MNI,1),1)*6];
ROI_name = {'Posterior cingulate cortex/Precuneus','PCC_DMN';'Medial prefrontal cortex','mPFC_DMN';'Left lateral parietal cortex','L_lPar_DMN';'Right lateral parietal cortex','R_lPar_DMN';'Left inferior temporal gyrus','L_IT_DMN';'Right inferior temporal gyrus','R_IT_DMN';'Medial dorsal thalamus','mdThal_DMN';'Left posterior cerebellum','L_pCERE_DMN';'Right posterior cerebellum','R_pCERE_DMN'};
files_conf = spm_select('FPListRec',main,'^sub-.*\_task-rest_bold_confounds.tsv$' );
for i=1:size(files_conf,1)
    fl = strcat(files_conf(i,:));
    [fpath,name,ext] = fileparts(fl);
    fnii = fullfile(fpath,[strrep(name,'confounds','space-MNI152NLin2009cAsym_preproc.nii')]);
    rsHRF_demo_UCLA_3D_ROI_deconv_sub_gamma3(fl,fnii,ROI,mainoutdir)
end

