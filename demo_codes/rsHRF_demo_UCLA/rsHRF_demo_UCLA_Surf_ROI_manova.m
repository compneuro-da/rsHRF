clc,clear
main = 'E:\data\ds000030_R1.0.5\derivatives\fmriprep';
sub= spm_load('E:\data\participants_ucla.tsv');
subid = cell2mat(sub.participant_id); 
subid = str2num(subid(:,5:end));
files_L = spm_select('FPListRec',main,'^sub-.*\_task-rest_bold_space-fsaverage5.L.func.gii$');
hrfdir = 'E:\data\HRF_outdir_surf_L_net7_Fourier';
isubk=1; 
for isub=1:size(files_L)%:-1:1
    fl = strcat(files_L(isub,:));
    [fpath,name,ext] = fileparts(fl);
    sub_id = str2num(fpath(end-9:end-5));
    idsub = find(sub_id == subid);
    if isempty(idsub)
        error(num2str(sub_id))
    end
    conf = strrep(fl,['space-fsaverage5.L.func.gii'],'confounds.tsv');

    confa = spm_load(conf);
    mfdr = confa.FramewiseDisplacement;
    mfdr(1) = 0;
    mfd0(isub,1) = mean(mfdr);
    if sub.ScannerSerialNumber(idsub)== 35343
        scan0  = 1;
    elseif sub.ScannerSerialNumber(idsub)== 35426
        scan0  = 2;
    end
    if mfd0(isub)<0.2  
        subject{isubk,1} = sub_id;
        mfd(isubk,1) = mfd0(isub);
        age(isubk,1) = sub.age(idsub);
        sex{isubk,1} = sub.gender{idsub};
        group{isubk,1} = sub.diagnosis{idsub};
        scans{isubk,1} = ['scanner' num2str(scan0)];
        Y{isubk,1} = sub.diagnosis{idsub};
        
        parastr = {'Height','FWHM','Time2peak'};
        fl0 = strrep(fl,'func\','');
        im = strrep(fl0,'.gii','_hrf.mat');
        im = spm_file(im,'prefix','T5to3AR1_');
        im = strrep(im,main,hrfdir);
        im = strrep(im,'\','/');
        hrfnii{isubk,1} = im;
        isubk=isubk+1;        
    end
end

if 1
    hrfdata= [];
    for i=size(hrfnii,1):-1:1
        load(hrfnii{i,1})
        hrfdata(i,:,:) = PARA;
    end
end

F = [];pvalue=[];
tic
Meas = table(parastr','VariableNames',{'HRFparameter'});

for i=1:size(hrfdata,3)
    dat = hrfdata(:,:,i);
    tall = table(Y,sex,scans,age,mfd,dat(:,1),dat(:,2),dat(:,3),'VariableNames',{'Group','Sex','Scanner','Age','meanFD','hrf1','hrf2','hrf3'});
    rm = fitrm(tall,'hrf1-hrf3~Group*Scanner+Sex+Age+meanFD','WithinDesign',Meas);

    ma = manova(rm);
    F(i,:) = ma.F;
    pvalue(i,:) = ma.pValue;
    
%     for j=1:3
%         tallh = table(Y,sex,scans,age,mfd,dat(:,j),'VariableNames',{'Group','Sex','Scanner','Age','meanFD','hrf'});
%         lm = fitlm(tallh,'hrf~Group*Scanner+Sex+Age+meanFD');
%         ah = anova(lm);
%         hrf_pvalue(i,j) = ah.pValue(1);
%     end
    mgs = multcompare(rm,'Group');
%     mgs = multcompare(rm,'Group','By','Scanner');
    posthoc_pvalue(i,:) = mgs.pValue;
    posthoc_Difference(i,:) = mgs.Difference;
end
fprintf('MANOVA Done, %5.2fs\n',toc);
fdr0 = mafdr(pvalue(:,8),'BHFDR',true);
rsn7={'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default'};
disp('Roy''s maximum root statistic:')
for i=1:7
    fprintf('%20s, p=%1.5f, p(FDR)=%1.5f\n',rsn7{i},pvalue(i,8),fdr0(i))
end
nnz(posthoc_pvalue<0.05)