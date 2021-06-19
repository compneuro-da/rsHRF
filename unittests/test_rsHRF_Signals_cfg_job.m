function tests = test_rsHRF_Signals_cfg_job
% Unit Tests for (Signal)ROI rsHRF estimation/deconvolution and FC analysis
tests = functiontests(localfunctions);

function test_Signals(testCase)
out_dir = tempdir;
conprefix = 'Conntest_';
deconvprefix = 'Deconvtest_';
filename = 'rsHRF_demo';
expected_output = {fullfile(out_dir,[deconvprefix,'combROI_',filename,'_job.mat'])
fullfile(out_dir,[deconvprefix,'combROI_',filename,'_hrf.mat'])
fullfile(out_dir,[conprefix,'combROI_',filename,'_deconv_Corr_Spearman.mat'])
};
fpath = fileparts(which('rsHRF.m'));
gii = fullfile(fpath,'unittests','rsHRF_demo.gii');
a = gifti(gii); 
dd = double(a.cdata');
sig1 = fullfile(out_dir,[filename,'.txt']);
sig2 = fullfile(out_dir,[filename,'.mat']);
save(sig1,'dd','-ascii');
save(sig2,'dd');
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Datasig(1).sigdata = {sig2};
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Datasig(1).name = 'dd';
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Datasig(2).sigdata = {sig1};
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Datasig(2).name = '';
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Denoising.generic = {};
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Denoising.Detrend = 0;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Denoising.BPF = {};
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.Denoising.Despiking = 0;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.hrfm = 2;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.TR = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.hrflen = 32;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.num_basis = NaN;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.mdelay = [4 8];
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.cvi = 0;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.fmri_t = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.fmri_t0 = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.thr = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.localK = 2;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.tmask = NaN;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.hrfdeconv = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.connectivity{1}.conn.data4conn = 2;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.connectivity{1}.conn.method = 6;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.connectivity{1}.conn.prefix = conprefix;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.outdir = {out_dir};
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.savedata.deconv_save = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.savedata.hrfmat_save = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.savedata.job_save = 1;
matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.prefix = deconvprefix;
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

job = load(expected_output{1}); job = job.job;
testCase.assertTrue(isfield(job,'prefix'));

job = load(expected_output{2});
testCase.assertEqual(size(job.PARA), [3 10]);

job = load(expected_output{3});
testCase.assertTrue(isfield(job,'M'));
testCase.assertEqual(size(job.M.Matrix_z), [10 10]);

try
    delete(fullfile(out_dir,[deconvprefix '*']))
    delete(fullfile(out_dir,[conprefix,'*']))
end

matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.connectivity = cell(1, 0);
for i=[1 3:7]
    matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.HRFE.hrfm = i;
    deconvprefix2 = ['Deconvtest_',num2str(i),'_'];
    matlabbatch{1}.spm.tools.rsHRF.sig_rsHRF.prefix = deconvprefix2;
    spm_jobman('run', matlabbatch);
    expected_output2 = fullfile(out_dir,[deconvprefix2,'combROI_',filename,'_hrf.mat']);
    job = load(expected_output2);
    testCase.assertEqual(size(job.PARA), [3 10]);
end
try
    delete(fullfile(out_dir,[deconvprefix,'*']))
end