function tests = test_rsHRF_Voxelwise_cfg_job
% Unit Tests for voxelwise rsHRF estimation/deconvolution and FC/GC analysis
tests = functiontests(localfunctions);

function test_voxelwise(testCase)
out_dir = tempdir;
conprefix = 'Conntest_';
deconvprefix = 'Deconvtest_';
filename = 'rsHRF_demo';
expected_output = {fullfile(out_dir,[deconvprefix,filename,'_job.mat'])
fullfile(out_dir,[deconvprefix,filename,'_hrf.mat'])
fullfile(out_dir,[conprefix,filename,'_deconv_pwGC.mat'])
fullfile(out_dir,[conprefix,'10_',filename,'_deconv_Olrm_Z_Pearson.nii'])
};
fpath = fileparts(which('rsHRF.m'));
fnii = fullfile(fpath,'unittests','rsHRF_demo.nii,1');
atlas_nii = fullfile(fpath,'unittests','rsHRF_demo_atlas.nii,1');
mask_nii= fullfile(fpath,'unittests','rsHRF_demo_mask.nii,1');
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
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.rmoutlier = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.mask = {''};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{1}.conn.data4conn = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{1}.conn.Seed_ROI = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{1}.conn.genericROI{1}.atlas = {atlas_nii};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{1}.conn.method = 4;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{1}.conn.prefix = conprefix;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.data4conn = 2;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.Seed_ROI = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.genericROI{1}.atlas = {atlas_nii};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.method = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.morder = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.ndinfo = [NaN NaN];
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity{2}.connG.prefix = conprefix;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.outdir = {out_dir};
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.deconv_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.hrfnii_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.prefix = deconvprefix;
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

job = load(expected_output{1}); job = job.job;
testCase.assertTrue(isfield(job,'prefix'));

job = load(expected_output{2});
testCase.assertEqual(size(job.PARA_rm), [3 27]);

job = load(expected_output{3});
testCase.assertTrue(isfield(job,'M'));
testCase.assertEqual(size(job.M.GC_Matrix), [10 10]);

testCase.assertTrue(exist(expected_output{4},'file') > 0);

try
    delete(fullfile(out_dir,[deconvprefix '*']))
    delete(fullfile(out_dir,[conprefix,'*']))
end


matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.rmoutlier = 0;
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.connectivity = cell(1, 0);
matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.mask = {mask_nii};
for i=[1 3:7]
    matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.HRFE.hrfm = i;
    deconvprefix2 = ['Deconvtest_',num2str(i),'_'];
    matlabbatch{1}.spm.tools.rsHRF.vox_rsHRF.prefix = deconvprefix2;
    spm_jobman('run', matlabbatch);
    expected_output2 = fullfile(out_dir,[deconvprefix2,filename,'_hrf.mat']);
    job = load(expected_output2);
    testCase.assertEqual(size(job.PARA), [3 1]);
end
try
    delete(fullfile(out_dir,[deconvprefix,'*']))
end
