function rsHRF_run_all_unittests()
% Execute all unit tests contained in ./rsHRF/unittests subfolder.
fpath = fileparts(which('rsHRF.m'));
cd(fullfile(fpath,'unittests'))
runtests