clc,clear
main = 'E:\data\ds000030_R1.0.5\derivatives\fmriprep';
mainoutdir = 'E:\data\HRF_outdir_WM';
niimask = 'E:\rsHRF_demo_UCLA\dependent_file\wm_mask095.nii';
files_conf = spm_select('FPListRec',main,'^sub-.*\_task-rest_bold_confounds.tsv$' );
for i=1:size(files_conf,1)
    fl = strcat(files_conf(i,:));
    [fpath,name,ext] = fileparts(fl);
    csf_mat = fullfile(fpath,'CSF_sig.mat');
    fnii = fullfile(fpath,[strrep(name,'confounds','space-MNI152NLin2009cAsym_preproc.nii')]);
    rsHRF_demo_UCLA_sub_gamma3_wm(fl,fnii,niimask,mainoutdir,csf_mat)
end