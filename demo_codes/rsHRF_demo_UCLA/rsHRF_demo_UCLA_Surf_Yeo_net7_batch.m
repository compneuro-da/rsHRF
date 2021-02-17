clc,clear
main = 'E:\data\ds000030_R1.0.5\derivatives\fmriprep'; % data dir
outdir = 'E:\data\HRF_outdir_surf_L_net7_Fourier'; % out dir

atlasgii = 'E:\rsHRF_demo_UCLA\dependent_file\lh.Yeo2011_7Networks_N1000.gii'; % atlas
files_conf = spm_select('FPListRec',main,'^sub-.*\_task-rest_bold_confounds.tsv$' );
for i=1:size(files_conf,1)
    fl = strcat(files_conf(i,:));
    [fpath,name,ext] = fileparts(fl);
    gnii = fullfile(fpath,[strrep(name,'confounds','space-fsaverage5.L.func.gii')]);
    rsHRF_demo_UCLA_sub_Fourier_surf_ROI(fl,gnii,outdir,atlasgii)
end