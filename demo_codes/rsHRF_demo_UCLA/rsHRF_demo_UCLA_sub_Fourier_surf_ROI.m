%-----------------------------------------------------------------------
% Job saved on 01-Oct-2020 17:30:46 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
function rsHRF_demo_UCLA_sub_Fourier_surf_ROI(fl,gnii,mainoutdir,atlasgii)
[fpath,name,~] = fileparts(fl);
id = strfind(name,'_task');
subid = name(1:id(1)-1);
a= spm_load(fl);
acomp =   [a.aCompCor00, a.aCompCor01, a.aCompCor02, a.aCompCor03, a.aCompCor04, a.aCompCor05];
Q1  = [a.X, a.Y, a.Z, a.RotX, a.RotY, a.RotZ];
HM = [Q1, [zeros(1,size(Q1,2));diff(Q1)]];
FD = a.FramewiseDisplacement; FD(1)=0;
if mean(FD)>0.2
    return
end
tmask = double(FD<0.3);
nui = [acomp HM];
txtfile = fullfile(fpath,[name,'.txt']);
save(txtfile,'nui','-ascii')
outdir = fullfile(mainoutdir,subid); 

matlabbatch={};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.images = {gnii};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.generic{1}.multi_reg = {txtfile};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.Detrend = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.BPF{1}.bands = [0.01 0.1];
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.Despiking = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.which1st = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.genericROI{1}.meshatlas = {atlasgii};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrfm = 5;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.TR = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrflen = 20;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.num_basis = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.mdelay = [4 8];
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.cvi = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.fmri_t = 5;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.fmri_t0 = 3;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.localK = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.tmask = tmask;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.mask = {''};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity = {};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.outdir = {outdir};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.deconv_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.prefix = 'T5to3AR1_';
spm('defaults', 'FMRI');
% spm_jobman('interactive', matlabbatch);
spm_jobman('run', matlabbatch);