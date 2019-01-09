clc,clear
rawcode_dir = fileparts(which('rsHRF_install_SPM.m'));
rsHRF_file = {
'wgr_get_parameters.m'  
'wgr_rshrf_estimation_canonhrf2dd_par2.m'
'wgr_rsHRF_FIR.m'
'wgr_rsHRF_global_para.m'
'rsHRF.m'
'rsHRF_install_SPM.m'
'rsHRF_toolbox.pptx'
'rsHRF.man'
'rsHRF_logo.png'
'spm_rsHRF.m'
'tbx_cfg_rsHRF.m'
'LICENSE'
'deleteoutliers.m'
'inpaint_nans3.m'
'knee_pt.m'
'README.md'
'demo_jobs.zip'
};
numfile = length(rsHRF_file);
toolbox_dir = fullfile(spm('Dir'),'toolbox','rsHRF');
try
    if all(toolbox_dir==rawcode_dir)
        disp('you had installed rsHRF.')
        return
    end
end
if ~exist(toolbox_dir,'dir')
    mkdir(toolbox_dir);
else
    try
        zip(fullfile(toolbox_dir,['rsHRF',date,'.zip']),rsHRF_file,toolbox_dir)
    end
end

for i=1:numfile
    file = fullfile(toolbox_dir,rsHRF_file{i});
    if exist(file,'file')
        delete(file)
    end
    copyfile(fullfile(rawcode_dir,rsHRF_file{i}),file)
end
unzip(fullfile(toolbox_dir,'demo_jobs.zip'),toolbox_dir)
disp('Done, rsHRF')
str='rsHRF';
disp(str)
eval(str)
% close
pause(3)
disp('Connectivity')
str2= 'rsHRF conn';
eval(str2)
pause(3)
close all