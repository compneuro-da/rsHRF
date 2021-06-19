function tests = test_rsHRF_SurfROI_cfg_job
% Unit Tests for (surface) ROI rsHRF estimation/deconvolution and FC/GC analysis
tests = functiontests(localfunctions);

function test_SurfROI(testCase)
out_dir = tempdir;
conprefix = 'Conntest_';
deconvprefix = 'Deconvtest_';
filename = 'rsHRF_demo';
expected_output = {fullfile(out_dir,[deconvprefix,filename,'_job.mat'])
fullfile(out_dir,[deconvprefix,filename,'_hrf.mat'])
fullfile(out_dir,[conprefix,filename,'_deconv_PCGC.mat'])
fullfile(out_dir,[conprefix,filename,'_deconv_Corr_PartialPearson.mat'])
};
fpath = fileparts(which('rsHRF.m'));
fgii = fullfile(fpath,'unittests','rsHRF_demo.gii');
atlas_gii = fullfile(fpath,'unittests','rsHRF_demo_atlas.gii');
mask_gii= fullfile(fpath,'unittests','rsHRF_demo_mask.gii'); 
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.images = {fgii};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.generic = {};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.Detrend = 0;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.BPF = {};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.Despiking = 0;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.Denoising.which1st = 3;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.genericROI{1}.meshatlas = {atlas_gii};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrfm = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.TR = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrflen = 32;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.num_basis = NaN;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.mdelay = [4 8];
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.cvi = 0;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.fmri_t = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.fmri_t0 = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.localK = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.tmask = NaN;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrfdeconv = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.mask = {''};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{1}.conn.data4conn = 3;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{1}.conn.method = 5;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{1}.conn.prefix = conprefix;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{2}.connG.data4conn = 2;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{2}.connG.method = 3;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{2}.connG.morder = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{2}.connG.ndinfo = [2 3];
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity{2}.connG.prefix = conprefix;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.outdir = {out_dir};
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.deconv_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.prefix = deconvprefix;
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

job = load(expected_output{1}); job = job.job;
testCase.assertTrue(isfield(job,'prefix'));

job = load(expected_output{2});
testCase.assertEqual(size(job.PARA), [3 5]);

job = load(expected_output{3});
testCase.assertTrue(isfield(job,'M'));
testCase.assertEqual(size(job.M.GC_Matrix), [5 5]);

testCase.assertTrue(exist(expected_output{4},'file') > 0);

try
    delete(fullfile(out_dir,[deconvprefix '*']))
    delete(fullfile(out_dir,[conprefix,'*']))
end


matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.connectivity = cell(1, 0);
matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.mask = {mask_gii};
for i=[1 3:7]
    matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.HRFE.hrfm = i;
    deconvprefix2 = ['Deconvtest_',num2str(i),'_'];
    matlabbatch{1}.spm.tools.rsHRF.SurfROI_rsHRF.prefix = deconvprefix2;
    spm_jobman('run', matlabbatch);
    expected_output2 = fullfile(out_dir,[deconvprefix2,filename,'_hrf.mat']);
    job = load(expected_output2);
    testCase.assertEqual(size(job.PARA), [3 5]);
end
try
    delete(fullfile(out_dir,[deconvprefix,'*']))
end

