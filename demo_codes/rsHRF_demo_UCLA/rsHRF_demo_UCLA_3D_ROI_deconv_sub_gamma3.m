function rsHRF_demo_UCLA_3D_ROI_deconv_sub_gamma3(fl,fnii,ROI,mainoutdir)
[fpath,name,~] = fileparts(fl);
id = strfind(name,'_task');
subid = name(1:id(1)-1);
gunzip([fnii,'.gz'])
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

matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.images = {[fnii,',1']};
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.Denoising.generic{1}.multi_reg = {txtfile};
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.Denoising.Detrend = 1;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.Denoising.BPF{1}.bands = [0.01 0.1];
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.Denoising.Despiking = 1;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.genericROI{1}.ROI = ROI;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.hrfm = 3;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.TR = 2;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.hrflen = 32;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.num_basis = 3;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.mdelay = [3 9];
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.cvi = 0;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.fmri_t = 5;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.fmri_t0 = 3;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.tmask = NaN;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.HRFE.tmask = tmask;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.mask = {''};
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.connectivity{1}.conn.data4conn = 3;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.connectivity{1}.conn.method = 5;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.connectivity{1}.conn.prefix = 'Conn_';
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.outdir = {outdir};
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.savedata.deconv_save = 1;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.savedata.job_save = 0;
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.prefix = 'Deconv_';
spm('defaults', 'FMRI');
spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);

delete(fnii)