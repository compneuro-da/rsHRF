clc,clear
main = 'E:\data\HRF_outdir_WM'; % data directory
mains  = 'E:\data\HRF_WM_HRFpara_3D'; % out directory
para = {'Height','Time2peak','FWHM'};
for k=1:3
    str = para{k};
    files = spm_select('FPListRec',main,['^T5to3AR1.*\_Olrm_',str,'.nii$']);
    for i=1:size(files,1)
        a = strcat(files(i,:));
        b = spm_file(a,'prefix','s');
        b = strrep(b,main,mains);
        fpath = fileparts(b);
        if ~exist(fpath)
            mkdir(fpath)
        end
        spm_smooth(a,b,[6 6 6])
        fprintf('%s\n',a)
    end
end