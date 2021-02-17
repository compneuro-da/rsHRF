function rsHRF_demo_UCLA_sub_gamma3_wm(fl,fnii,niimask,mainoutdir,csf_mat)
[fpath,name,~] = fileparts(fl);
id = strfind(name,'_task');
subid = name(1:id(1)-1);
gunzip([fnii,'.gz'])
a= spm_load(fl);
b = load(csf_mat);
acomp =   [b.csf_acompcor];
Q1  = [a.X, a.Y, a.Z, a.RotX, a.RotY, a.RotZ];
HM = [Q1, [zeros(1,size(Q1,2));diff(Q1)]];
FD = a.FramewiseDisplacement; FD(1)=0;
tmask = double(FD<0.3);
if mean(FD)>0.2
    return
end
nui = [acomp HM];
txtfile = fullfile(fpath,[name,'_csf.txt']);
save(txtfile,'nui','-ascii')
outdir = fullfile(mainoutdir,subid); 
rsHRF_gamma3HRF_batch_job([fnii,',1'],txtfile,tmask,niimask,outdir)
delete(fnii)

function rsHRF_gamma3HRF_batch_job(fnii,txtfile,tmask,nii_mask,outdir)
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.images = {fnii};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.generic{1}.multi_reg = {txtfile};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.Detrend = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.BPF{1}.bands = [0.01 0.1];
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.Despiking = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrfm = 3;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.TR = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrflen = 32;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.mdelay = [4 8];
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.cvi = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.num_basis = 3;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.fmri_t = 5;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.fmri_t0 = 3;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.localK = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.tmask = tmask;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.rmoutlier = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.mask = {nii_mask};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity = {};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.outdir = {outdir};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.deconv_save = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfnii_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.prefix = 'T5to3AR1_';
spm('defaults', 'FMRI');
spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);