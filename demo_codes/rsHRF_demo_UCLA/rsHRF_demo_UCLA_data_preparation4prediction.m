clc,clear,close all
sub= spm_load('E:\data\participants_ucla.tsv');
subid = cell2mat(sub.participant_id); subid = str2num(subid(:,5:end));
main = 'E:\data\ds000030_R1.0.5\derivatives\fmriprep';
hrfdir = 'E:\data\UCLA_ROI9_FC';
files_conf = spm_select('FPListRec',main,'^sub-.*\_task-rest_bold_confounds.tsv$' );

for scanner=1:2
    age=[]; gender=[]; scan=[]; mfd =[]; datz =[]; datzd=[];
    isubk=1; 
    for isub=1:size(files_conf,1)
        fl0 = strcat(files_conf(isub,:));
        [~,name,~] = fileparts(fl0);
        sub_id = str2num(name(5:9));
        fpath = fileparts(strrep(fileparts(fl0),main,hrfdir)); 

        idsub = find(subid == sub_id);
        fprintf('%d,%d\n',isubk,idsub)

        if sub.gender{idsub}=='F'
            gender0  = 1;
        elseif sub.gender{idsub}=='M'
            gender0  = 2;
        else
            error('unknown gender information')
        end

        if sub.ScannerSerialNumber(idsub)== 35343
            scan0  = 1;
        elseif sub.ScannerSerialNumber(idsub)== 35426
            scan0  = 2;
        end

        confa = spm_load(fl0);
        fd0 = confa.FramewiseDisplacement;
        fd0(1) = 0;
        fd = mean(fd0);
        
        if fd<0.2 & sub.diagnosis{idsub}(1)=='C'  & scan0==scanner
            age(isubk,1) = sub.age(idsub);
            gender(isubk,1) = gender0;
            scan(isubk,1) = scan0;
            Y{isubk,1} = sub.diagnosis{idsub};
            mfd(isubk,1) = fd;
            f0 = spm_select('FPListRec',fpath,'^Conn.*\_preproc_Corr_PartialPearson.mat$');
            load(f0)
            if ~exist('idm','var')
                idm = find(tril(ones(size(M.Matrix_z)),-1));  
            end
            datz(isubk,:) = M.Matrix_z(idm);
            clear M
            f1 = spm_select('FPListRec',fpath,'^Conn.*\_preproc_deconv_Corr_PartialPearson.mat$' );
            load(f1)
            datzd(isubk,:) = M.Matrix_z(idm);
            clear M
            isubk=isubk+1;
        end
    end
    
    cova = [gender mfd];
    all_matr  = datz;
    all_matd  = datzd;
    save(['FC',num2str(length(age)),'.mat'],'cova','all_*','age')
end