function tests = test_rsHRF_Viewer_cfg_job
% Unit Tests for voxelwise rsHRF display
tests = functiontests(localfunctions);

function test_Viewer(testCase)
out_dir = tempdir;
deconvprefix = 'Deconvtest_';
filename = 'rsHRF_demo';
expected_output = {fullfile(out_dir,[deconvprefix,filename,'_job.mat'])
fullfile(out_dir,[deconvprefix,filename,'_hrf.mat'])
};
fpath = fileparts(which('rsHRF.m'));
fnii = fullfile(fpath,'unittests','rsHRF_demo.nii,1');
atlas_nii = fullfile(fpath,'unittests','rsHRF_demo_atlas.nii,1');
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.images = {fnii};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.generic = {};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.Detrend = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.BPF = {};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.Despiking = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.Denoising.which1st = 3;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrfm = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.TR = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrflen = 32;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.num_basis = 3;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.mdelay = [4 8];
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.cvi = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.fmri_t = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.fmri_t0 = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.localK = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.tmask = NaN;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrfdeconv = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.rmoutlier = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.mask = {''};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.outdir = {out_dir};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.deconv_save = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfnii_save = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.prefix = deconvprefix;
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

matlabbatch={};
matlabbatch{1}.spm.tools.rsHRF.display_HRF.underlay_nii = {atlas_nii};
matlabbatch{1}.spm.tools.rsHRF.display_HRF.stat_nii = {atlas_nii};
matlabbatch{1}.spm.tools.rsHRF.display_HRF.HRF_mat = {expected_output([2 2 2])};
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);
try
    delete(fullfile(out_dir,[deconvprefix '*']))
end
close all