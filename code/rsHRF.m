function rsHRF(varargin)
%% resting state BOLD-fMRI HRF deconvolution and connectivity analysis
if nargin==2
    job = varargin{1};  flag = varargin{2};  
    if ~isfield(job,'para_global')
        para_global = wgr_rsHRF_global_para;
        job.para_global = para_global; clear para_global
    end
    switch flag
        case 'display'
            rsHRF_viewer(job);
        case 'volume'
            % delete temporary files (generated wm/csf/brainmask) 
            flag_delete = job.para_global.delete_files; 
            
            %input images, cell file
            nii_file = job.images; 
            
            %get 1st file name (out name)
            ref_nii = nii_file{1}; %funtional space template
            [fpath,raw_outname,~] = fileparts(ref_nii);
            job.ref_nii = ref_nii;
            job.raw_outname = raw_outname;
            
            %out directory
            outdir = job.outdir{1};
            if isempty(outdir)
                outdir = fpath;
            else
                if ~exist(outdir,'dir')
                    mkdir(outdir)
                end
            end
            
            % NIfTI data 
            [data,v] = wgr_read_NIfTI_job(job);

            % explicit mask
            nii_mask = job.mask{1};
            if isempty(nii_mask)
                nii_mask = job.para_global.mask_nii ;
            end      
            if ~isempty(nii_mask)  %brain mask
                nii_masknew  = spm_file(nii_mask,'prefix',['bin_',raw_outname]);
                spm_imcalc({ref_nii;nii_mask},nii_masknew,'(i1~=0).*(i2~=0)',{0,0,0,16});
                v0 = spm_vol(nii_masknew);
                datam = spm_read_vols(v0);
                datam(isnan(datam))=0;
                voxel_ind = find(datam~=0); clear datam
                if flag_delete
                    delete(v0.fname)
                end
            else % var>0 as mask
                data_var = nanvar(data,0,1);
                voxel_ind = find(data_var>0); clear data_var
            end
            fprintf('# %d voxels inside defaut mask ...\n',length(voxel_ind))

           %% denoise
            data = wgr_Denoising_job(job,data);
            
            % check ROI information
            if  isfield(job,'genericROI')    
                ROI = job.genericROI;
            else
                ROI = [];
            end
            
            if ~isempty(ROI)
                flag_ROI = 1;
                [mat,atlas,nii]= wgr_check_ROI(ROI, ref_nii,flag_delete);                
                roinfo = [mat;atlas;nii];
                nvar = size(roinfo,1);
                Nscans = size(data,1);
                data_ROI = zeros(Nscans,nvar);
                for j=1:nvar
                    data_ROI(:,j) = mean(data(:,roinfo{j,1}),2);
                end
            else
                flag_ROI=0;
            end
            
           %% data for further analysis
            if flag_ROI
                data = data_ROI;
                
                if isfield(job,'prefix')
                    if job.savedata.job_save
                        save(fullfile(outdir,[job.prefix,raw_outname, '.mat']), 'roinfo','-v7.3');
                    end
                else
                    if job.job_save
                        save(fullfile(outdir,[raw_outname, '_roinfo.mat']), 'roinfo');
                    end
                end
            else
                data = data(:,voxel_ind);
            end
           
            flag_nii_gii = 1; %NIfTI
            if isfield(job,'HRFE') % deconvolution
                wgr_deconv_job(job,data,flag_ROI, outdir, v, voxel_ind, flag_nii_gii);
            else
                wgr_conn_job(job,data,flag_ROI, outdir, v, voxel_ind, flag_nii_gii)
            end
            
        case 'mesh'    
            %input images, cell file
            gii_file = job.images; 
            [fpath,raw_outname,~] = fileparts(gii_file{1});
            job.raw_outname = raw_outname;
            job.rmoutlier = 0;
            %out directory
            outdir = job.outdir{1};
            if isempty(outdir)
                outdir = fpath;
            else
                if ~exist(outdir,'dir')
                    mkdir(outdir)
                end
            end
            % GIfTI data 
            data = wgr_read_GIfTI_job(job);
            v.dim = [size(data,2), 1];
            % explicit mask
            gii_mask = job.mask{1};

            if ~isempty(gii_mask)  %brain mask
                datam = gifti(gii_mask);
                datam.cdata(isnan(datam.cdata))=0;
                vertex_ind = find(datam.cdata~=0); clear datam
            else % var>0 as mask
                data_var = nanvar(data,0,1);
                vertex_ind = find(data_var>0); clear data_var
            end
            fprintf('# %d vertexs inside the mask ...\n',length(vertex_ind))
            %% Denoise
            data = wgr_Denoising_job(job,data);

            % check ROI information
            if  isfield(job,'genericROI')    
                ROI = job.genericROI;
            else
                ROI = [];
            end

            if ~isempty(ROI)
                flag_ROI = 1;
                [~,~,~,atlasmeh,gii]= wgr_check_ROI(ROI);                
                roinfo = [atlasmeh; gii];
                nvar = size(roinfo,1);
                Nscans = size(data,1);
                data_ROI = zeros(Nscans,nvar);
                for j=1:nvar
                    data_ROI(:,j) = nanmean(data(:,roinfo{j,1}),2);
                end
            else
                flag_ROI=0;
            end

            %% data for further analysis            
            if flag_ROI
                data = data_ROI;

                if isfield(job,'prefix')
                    if job.savedata.job_save
                        save(fullfile(outdir,[job.prefix,raw_outname, '.mat']), 'roinfo','-v7.3');
                    end
                else
                    if job.job_save
                        save(fullfile(outdir,[raw_outname, '_roinfo.mat']), 'roinfo');
                    end
                end
            else
                data = data(:,vertex_ind);
            end
            
            flag_nii_gii = 2; %GIfTI
            if isfield(job,'HRFE') % deconvolution
                wgr_deconv_job(job,data,flag_ROI, outdir, v, vertex_ind, flag_nii_gii);
            else
                wgr_conn_job(job,data,flag_ROI, outdir, v, vertex_ind, flag_nii_gii)
            end
        case 'sig'
            wgr_ROI_sig_job(job)         
    end
    return
else
    code_path = fileparts(which('rsHRF.m')) ;
    a = imread(fullfile(code_path,'rsHRF_logo.png'));
    pos = [0.4   0.3    0.2    0.3236];
    if ismac       % Code to run on Mac plaform
    elseif isunix  % Code to run on Linux plaform
        pos = [0.4   0.3    0.23    0.35];
    elseif ispc    % Code to run on Windows platform
    end

    S.fig = figure('Visible','on',...
        'numbertitle','off',...    
        'menubar','none',...         
        'units','normalized',...
        'color','w',...
        'position',pos,...[563    98   480   360],... %0.4   0.3    0.15    0.38
        'name',['rsHRF v2.2(',getenv('USERNAME'),')'],...
        'resize','off');
    
    axes('parent',S.fig,'units','normalized','position',[0.15 0.45 [1 0.7825]*0.76]); %[0 0.59 1 0.43]
    imagesc(a, 'AlphaData', 0.6);
    axis off; 
    
    %% ROI signal HRF deconvolution
    S.pb(1) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.05,0.05,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','Signals',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.sig_rsHRF'');',...
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');

    %% ROI-wise HRF deconvolution
    S.pb(2) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.05,0.2,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','ROIs-volume',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.ROI_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');
    
    %% (Surface) ROI-wise HRF deconvolution
    S.pb(3) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.55,0.2,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','ROIs-surface',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.SurfROI_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');
    
    %% voxelwise HRF deconvolution
    S.pb(4) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.05,0.35,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','Voxels',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.vox_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');
    %% vertexwise HRF deconvolution
    S.pb(5) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.55,0.35,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','Vertices',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.mesh_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');
    
    S.pb(6) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.55,0.05,0.4,0.12],...
        'style','pushbutton',...'pushbutton',...
        'string','Display',...        
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.display_HRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.4,... 
        'fontweight','bold');
    
    if nargin==1
        Modality=upper(varargin{1});
        switch Modality
            case 'CONN'  
                set(S.pb(1),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.sig_conn'');'); 
                set(S.pb(2),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.ROI_conn'');');
                set(S.pb(3),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.SurfROI_conn'');');
                set(S.pb(4),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.vox_conn'');');
                set(S.pb(5),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.mesh_conn'');');
                set(S.pb(6),'enable','off')
                set(S.pb(:),'foregroundcolor',0*[1 1 1]) %0.3 0.6 0
                text(80,360,'Connectivity Analysis','Color',[0.3 0.6 0],'Fontsize',13,'Fontweight','bold');
            otherwise
                text(80,360,'Resting State HRF','Color',[1 0.6 0],'Fontsize',13,'Fontweight','bold');
        end
    else
        text(80,360,'Resting State HRF','Color',[1 0.6 0],'Fontsize',13,'Fontweight','bold');
    end
                                                                                                                                                      
end

function [data,v] = wgr_read_NIfTI_job(job)
%read NIfTI data 
nii_file = job.images; 
vol_th = job.para_global.volume_threshold;
data = [];
fprintf('Reading NIfTI Data ...\n')
Nscans = length(nii_file);
if Nscans<vol_th
    if Nscans==1
        [pth0,nam0,ext0] = spm_fileparts(nii_file{1});
        tmp = fullfile(pth0,[nam0,ext0]);
        v  = spm_vol(tmp);
        Nscans = length(v);
        if Nscans>vol_th
            for i=1:Nscans
                nii_file{i,1}  = [v(1).fname,',',num2str(i)];
            end
        end
    else
        error(['Please check your input data: ',num2str(Nscans),' volumes ? (set frames "inf" for 4D NIfTI data)'])
    end
end
for i=Nscans:-1:1
    v(i) = spm_vol(nii_file{i});
    tmp  = spm_read_vols(v(i));
    data(i,:) = tmp(:);
end       

function data = wgr_read_GIfTI_job(job)
%read NIfTI data 
gii_file = job.images; 
vol_th = job.para_global.volume_threshold;
data = [];
fprintf('Reading GIfTI Data ...\n')
Nscans = length(gii_file);
if Nscans<vol_th
    if Nscans==1
        [pth0,nam0,ext0] = spm_fileparts(gii_file{1});
        tmp = fullfile(pth0,[nam0,ext0]);
        v  = gifti(tmp);
        data = v.cdata'; %be careful
        [nobs,nvar] = size(data);
        fprintf('(Check) input data: %d time points, #%d vertices\n',nobs,nvar) 
    else
        error(['More than one file, please check your input data: ',num2str(Nscans)]) 
    end
end

function [mat,atlas,nii,atlasmesh,gii]= wgr_check_ROI(ROI, funii, flag_delete)
mat = {};  atlas = {}; nii = {}; atlasmesh = {}; gii = {};
k1 =1 ; k2 = 1; k3 = 1; k4 = 1; k5 = 1;
for i=1:length(ROI)
    tmp  = ROI{i};
    flag = isfield(tmp,'ROI');
    if flag
        for j=1: size(tmp.ROI,1)
            tmp2 = tmp.ROI(j,:);
            xY = [];
            xY.def='sphere';
            xY.xyz=tmp2(1:3)';
            xY.spec=tmp2(4);
            [xY, XYZmm, id] = spm_ROI(xY, funii);
            if isempty(id)
                error(sprintf('\n\n====(ROI)Error Information====\nno voxel found in sphere with center (%3.1f, %3.1f, %3.1f) and radius %3.1f.\nTry to increase the radius???\n',tmp2));
            else
                fprintf('#%d voxels found in sphere with center (%3.1f, %3.1f, %3.1f) and radius %3.1f.\n',length(id),tmp2)
            end
            mat{k1,1} = id;
            mat{k1,2} = tmp2;
            mat{k1,3} = '';
            k1 = k1+1;
            continue;
        end
    end
        
    flag = isfield(tmp,'atlas');
    if flag
        for j=1: length(tmp.atlas)
            tmp2 = strcat(tmp.atlas{j});
            if ~isempty(tmp2)
                atlas{k2,1} = tmp2;
                niinew  = spm_file(tmp2,'prefix','bin_');
                spm_imcalc({funii;tmp2},niinew,'(i1~=0).*i2',{0,0,0,16});
                v = spm_vol(niinew);
                datam = spm_read_vols(v); datam(isnan(datam)) = 0;
                idu = unique(datam(:)); idu(idu==0)=[];
                for k=1:length(idu)
                    atlas{k2,1} = find(datam==idu(k)); 
                    atlas{k2,2} = tmp2; 
                    atlas{k2,3} = idu(k); 
                    k2 = k2+1;
                end            
                clear datam
                if flag_delete
                    delete(v.fname)
                end
            end
        end
        continue;
    end
    
    flag = isfield(tmp,'images');
    if flag
        for j=1: length(tmp.images)
            tmp2 = strcat(tmp.images{j});      
            if ~isempty(tmp2)
                nii{k3,1} = tmp2;
                niinew  = spm_file(tmp2,'prefix','bin_');
                spm_imcalc({funii;tmp2},niinew,'(i1~=0).*i2',{0,0,0,16});
                v  =spm_vol(niinew);
                datam = spm_read_vols(v);
                nii{k3,1} = find(datam); clear datam
                nii{k3,2} = tmp2; 
                nii{k3,3} = ''; 
                if flag_delete
                    delete(v.fname)
                end
                k3 = k3+1;
            end
        end
        continue;
    end
    
    flag = isfield(tmp,'meshmask');
    if flag
        for j=1: length(tmp.meshmask)
            tmp2 = strcat(tmp.meshmask{j});      
            if ~isempty(tmp2)
                gii{k4,1} = tmp2;
                datam = gifti(tmp2);
                gii{k4,1} = find(datam.cdata); clear datam
                gii{k4,2} = tmp2; 
                gii{k4,3} = ''; 
                k4 = k4+1;
            end
        end
        continue;
    end
    
    flag = isfield(tmp,'meshatlas');
    if flag
        for j=1: length(tmp.meshatlas)
            tmp2 = strcat(tmp.meshatlas{j});
            if ~isempty(tmp2)
                atlasmesh{k5,1} = tmp2;
                datam = gifti(tmp2);
                idu = unique(datam.cdata(:)); idu(idu==0)=[]; 
                for k=1:length(idu)
                    atlasmesh{k5,1} = find(datam.cdata==idu(k)); 
                    atlasmesh{k5,2} = tmp2; 
                    atlasmesh{k5,3} = idu(k); 
                    k5 = k5+1;
                end            
                clear datam
            end
        end
        continue;
    end
end

function [data_txt,var_name]= wgr_check_ROIsig(ROI)
nROI = length(ROI);
data_txt = cell(nROI,2);
var_name  = cell(nROI,1);
for i=1:nROI
    tmp2 = ROI(i).sigdata{1};
    [~,~,ext] = fileparts(tmp2);
    if strfind(ext,'txt')
        data_txt{i,1} = load(tmp2);
        data_txt{i,2} = tmp2;
    elseif strfind(ext,'mat')
        tmp3 = load(tmp2);
        var_name{i} = strcat(ROI(i).name);
        try
            eval(['data_txt{i,1} = tmp3.',var_name{i},';']);
        catch
            tmp4 = fieldnames(tmp3);
            if length(tmp4)==1
                eval(['data_txt{i,1} = tmp3.',tmp4{1},';']);
            else
                error(['there are more than one variable in mat file, and no variable ''',var_name,'" in mat file']);
            end
        end
        data_txt{i,2} = tmp2;
    end
end

function data = wgr_Denoising_job(job,data)
Nscans = size(data,1);
flag_delete = job.para_global.delete_files; % delete temporary files (generated wm/csf/brainmask) 
covariates = job.Denoising.generic;
if ~isempty(covariates)
    fprintf('Reading Covariates ...\n')
    [txt,mat,nii]= wgr_check_covariates(covariates); 
    if nargin==3 % NIFTI
        if ~isempty(nii) 
            nii2 = cell(1,length(nii)); 
            for i=1:length(nii)
                nii2{1,i} = spm_file(nii{i,1},'prefix',['bin_',job.raw_outname]);
                spm_imcalc({job.ref_nii;nii{i,1}},nii2{1,i},'(i1~=0).*(i2>0.9)',{0,0,0,16});
                v0 = spm_vol(nii2{1,i});
                datam = spm_read_vols(v0);
                if flag_delete
                    delete(v0.fname)
                end
                id = find(datam==1); clear datam
                if isempty(id)
                    error(['no survived voxels (>0.9) in ',nii{i,1}])
                else
                    if nii{i,2}                                
                        numcomps = nii{i,2};
                        z = data(:,id);
                        z(isnan(z))=0;
                        m = size(z,2) ;
                        if m<numcomps
                            error(['# components = ',num2str(m),', < ', num2str(numcomps)])
                        else
                            pc = pca(z);
                        end
                        nii2{1,i} = z*pc(:,1:numcomps);
                    else
                        nii2{1,i} = nanmean(data(:,id),2);
                    end
                end
            end
            nii = cell2mat(nii2);
            clear nii2 datam
        else
            nii = [];
        end
        nuisance = [txt mat nii];
    elseif nargin==2
        nuisance = [txt mat];
    end
else
    nuisance = [];
end

%% Denoising parameters
p_detrend = job.Denoising.Detrend; 
if p_detrend
    nuisance_detrend = zeros(Nscans,p_detrend);
    for i = 1:p_detrend
        d = (1:Nscans).^i;
        nuisance_detrend(:,i+1) = d(:)./(Nscans.^i);
    end
else
    nuisance_detrend= [];
end

if ~isempty(nuisance) && ~isempty(nuisance_detrend)
    fprintf('Removing nuisance regressors & Detrend................. \n ')
elseif isempty(nuisance) && ~isempty(nuisance_detrend)
    fprintf('Detrend................. \n ')
elseif  ~isempty(nuisance) &&  isempty(nuisance_detrend)
    fprintf('Removing nuisance regressors................. \n ')
end
nuisance = [nuisance_detrend nuisance];
if ~isempty(nuisance)
    nuisance = bsxfun(@minus, nuisance, sum(nuisance, 1)./Nscans);% Remove the mean
    nuisance = [orth(nuisance) ones(Nscans, 1)];
    mu   = sum(data, 1)./Nscans;    % mean 
    data = data - nuisance*(nuisance\data);
    data = bsxfun(@plus, data, mu);
end

BPF = job.Denoising.BPF;
[Bands,TR] = wgr_check_bpf(BPF);
if ~isempty(Bands)
    fprintf('Frequency filter.................... \n ')
    if isnan(TR)
        if isfield(job,'HRFE')
            TR = job.HRFE.TR ; 
        else
            error('Please add TR information for Filtering')
        end
    else
        if isfield(job,'HRFE')
            if TR ~= job.HRFE.TR ; 
                error('Different TR for HRF estimation and Filtering')
            end
        end
    end       
    data = wgr_band_filter(data,TR,Bands);
end

if job.Denoising.Despiking
    fprintf(' Temporal despiking with a hyperbolic tangent squashing function... \n')
    data = wgr_despiking(data);
end

function [txt,mat,nii]= wgr_check_covariates(covariates)
txt={};  mat = {};  nii  = {};
k1 =1 ; k2 = 1; k3 = 1;
for i=1:length(covariates)
    tmp  = covariates{i};
    flag = isfield(tmp,'multi_reg');
    if flag
        for j=1: length(tmp.multi_reg)
            tmp2 = strcat(tmp.multi_reg{j});
            [~,~,ext] = fileparts(tmp2);
            if ~isempty(strfind(ext,'txt'))
                txt{1,k1} = load(tmp2);
                k1 = k1+1;
            elseif ~isempty(strfind(ext,'mat'))                
                tmp3  = load(tmp2);
                try
                    mat{1,k2} = tmp3.R;
                catch
                    tmp4 = fieldnames(tmp3);
                    if length(tmp4)==1
                        eval(['mat{1,k2} = tmp3.',tmp4{1},';']);
                    else
                        error('there are more than one variable in your mat file');
                    end
                end
                k2 = k2+1;
            else
                error('Unknown File Types')
            end
            continue;
        end
    end
        
    flag = isfield(tmp,'imcov');
    if flag
        for j=1: length(tmp.imcov.images)
            tmp2 = strcat(tmp.imcov.images{j});
            nii{k3,1} = tmp2;
            if length(tmp.imcov.mpc)==length(tmp.imcov.images)
                nii{k3,2} = tmp.imcov.mpc(j);
            elseif length(tmp.imcov.mpc)==1
                nii{k3,2} = tmp.imcov.mpc(1);
            end
            k3 = k3+1;
        end
        continue;
    end
end
txt = cell2mat(txt);
mat = cell2mat(mat);

function x1 = wgr_despiking(x1)
my=repmat(median(x1,1),[size(x1,1),1]);
sy=repmat(4*median(abs(x1-my)),[size(x1,1),1]);
x1=my+sy.*tanh((x1-my)./max(eps,sy));

function [Bands,TR] = wgr_check_bpf(BPF)
if length(BPF)==1
    try
        Bands = BPF{1}.bands;
        TR = nan;
    catch
        error('please check your setting for band-pass filtering \n')
    end
elseif length(BPF)==2
    try
        Bands = BPF{1}.bands;
        TR = BPF{2}.TR;
    catch
        try
            Bands = BPF{2}.bands;
            TR = BPF{1}.TR;
        catch
            error('please check your setting for band-pass filtering \n')
        end
    end
else
    Bands =[]; TR=nan;
end

function  x = wgr_band_filter(x, TR,Bands,m)

if nargin<4
    m = 5000; %block size
end
nvar = size(x,2);
nbin = ceil(nvar/m);
for i=1:nbin
    if i~=nbin
        ind_X = (i-1)*m+1:i*m ;   
    else
        ind_X = (i-1)*m+1:nvar ;
    end
    x1 = x(:,ind_X);
    x1 = conn_filter(TR,Bands,x1,'full') +...%mean-center by default
     repmat(mean(x1,1),size(x1,1),1);    %mean back in
    x(:,ind_X) = x1;
end
clear x1;

function [y,fy]=conn_filter(rt,filter,x,option)
% from conn toolbox
USEDCT=true;
if nargin<4, option='full'; end
if strcmpi(option,'base'), Nx=x; x=eye(Nx); end
if USEDCT % discrete cosine basis
    Nx=size(x,1);
    fy=fft(cat(1,x,flipud(x)),[],1);
    f=(0:size(fy,1)-1);
    f=min(f,size(fy,1)-f);
    switch(lower(option))
        case {'full','base'}
            idx=find(f<filter(1)*(rt*size(fy,1))|f>=filter(2)*(rt*size(fy,1)));
            %idx=idx(idx>1);
            fy(idx,:)=0;
            k=1; %2*size(fy,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(fy,1)-numel(idx));
            y=real(ifft(fy,[],1))*k;
            y=y(1:Nx,:);
        case 'partial'
            idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
            %if ~any(idx==1), idx=[1,idx]; end
            y=fy(idx,:);
    end
else % discrete fourier basis
    fy=fft(x,[],1);
    f=(0:size(fy,1)-1);
    f=min(f,size(fy,1)-f);
    switch(lower(option))
        case 'full'
            idx=find(f<filter(1)*(rt*size(fy,1))|f>=filter(2)*(rt*size(fy,1)));
            %idx=idx(idx>1);
            fy(idx,:)=0;
            k=1; %2*size(fy,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(fy,1)-numel(idx));
            y=real(ifft(fy,[],1))*k;
        case 'partial'
            idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
            %if ~any(idx==1), idx=[1,idx]; end
            y=fy(idx,:);
    end
end

function out = wgr_deconv_job(job,data,flag_ROI, outdir, v, smask_ind, flag_nii_gii)
[Nscans,nvar_data] = size(data);
name = job.raw_outname;
flag_pval_pwgc = job.para_global.flag_pval_pwgc;
flag_save_psc = job.para_global.flag_save_psc;

%% HRF deconvolution
para.TR = job.HRFE.TR ; 
%%% the following parameter (upsample grid) can be > 1 only for Canonical/Fourier/Gamma.
%%% Set = 1 for FIR/sFIR
para.T  = job.HRFE.fmri_t; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid
para.T0 = job.HRFE.fmri_t; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
if para.T==1
    para.T0 = 1;
end
min_onset_search = min(job.HRFE.mdelay); % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = max(job.HRFE.mdelay); % maximum delay allowed between event and HRF onset (seconds)
para.dt  = para.TR/para.T; % fine scale time resolution.
para.order = job.HRFE.num_basis;
para.AR_lag = job.HRFE.cvi; % AR(1) noise autocorrelation.
para.thr = job.HRFE.thr; % (mean+) para.thr*standard deviation threshold to detect event.
para.len = job.HRFE.hrflen; % length of HRF, in seconds
para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
para.localK = job.para_global.localK;
if isnan(job.HRFE.tmask(1))
    temporal_mask = true(Nscans,1); % temporal mask.
else
    temporal_mask = logical(job.HRFE.tmask);
end


if job.HRFE.hrfm<6
    hrfmn = {
    'Canonical HRF (with time derivative)'
    'Canonical HRF (with time and dispersion derivatives)'              
    'Gamma functions'
    'Fourier set'
    'Fourier set (Hanning)'}; % 
    para.name = hrfmn{job.HRFE.hrfm};
    tic
    [beta_hrf, bf, event_bold] = rsHRF_estimation_temporal_basis(data,para, temporal_mask);
    hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
    hrf_baseline = beta_hrf(1+size(bf,2),:); %HRF baseline value for PSC calculation.
else 
    tic
    if  job.HRFE.hrfm==6
        para.estimation = 'FIR';
    elseif job.HRFE.hrfm==7
        para.estimation = 'sFIR';
    end
    para.T=1; % this needs to be = 1 for FIR
    [beta_hrf,event_bold] = rsHRF_estimation_FIR(data,para, temporal_mask);
    hrfa = beta_hrf(1:end-2,:); %HRF
    hrf_baseline = beta_hrf(end-1,:); %HRF baseline value for PSC calculation.
end
 
nvar = size(hrfa,2);
PARA = zeros(3,nvar);
event_number = nan(nvar,1);
parfor ii=1:nvar
    hrf1 = hrfa(:,ii);
    PARA(:,ii) = wgr_get_parameters(hrf1,para.TR/para.T);% estimate HRF parameter
    event_number(ii)=length(event_bold{1,ii});
end
fprintf('\nDone HRF estimation %8.2f seconds\n',toc)
 
if job.savedata.hrfmat_save
    fprintf('Save HRF parameters....\n');
    %% save HRF parameters
    save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'para', 'beta_hrf', 'hrfa', 'event_bold', 'event_number','PARA','-v7.3');
    if exist('bf','var')
        save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'bf','-append')
    end
end
 
 if ~isempty(job.connectivity)
     [connroinfo,conndata_flag] = wgr_conn_check(job);
 else
     connroinfo={}; conndata_flag=0;
 end
    
%% HRF deconvolution -->signals
%whether save the deconvolved signal or perform analysis on deconvolved data
if job.savedata.deconv_save  || any(conndata_flag~=1) 
    fprintf('Deconvolving HRF ...\n');
    tic
%     T = round(para.len/TR);
    if para.T>1
        hrfa_TR = resample(hrfa,1,para.T);
    else
        hrfa_TR = hrfa;
    end
    data_deconv = nan(size(data));
    parfor ii=1:nvar
        hrf=hrfa_TR(:,ii);
        H=fft([hrf; zeros(Nscans-length(hrf),1)]);
        M=fft(data(:,ii));
        data_deconv(:,ii) = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));
    end
    data_deconv = real(data_deconv); %...
    toc
end
 
if ~flag_ROI % 3D volume/ 2D surface data
    if flag_nii_gii==1
        v0 = v(1);   v0.dt = [16,0]; v0.n = [1 1];
        ext_nii_gii = '.nii';
    else
        v0 = v;
        ext_nii_gii = '.gii';
    end
    HRF_para_str = {'Height', 'Time2peak','FWHM'};
    dat= zeros(v0.dim);
    if job.savedata.hrfnii_save
        %write nifti maps of the three HRF parameters
        for i=1:3
            fname = fullfile(outdir,[job.prefix,name,'_',HRF_para_str{i},ext_nii_gii]);
            eval(['out.',HRF_para_str{i},'{1} =','''',fname,''';'])
            dat(smask_ind)=PARA(i,:);
            wgr_write_file(fname,dat,flag_nii_gii,v0)
        end
        
        % write height PSC
        fname = fullfile(outdir,[job.prefix,name,'_Height_PSC',ext_nii_gii]);
        out.Height_PSC{1} = fname;
        dat(smask_ind)=PARA(1,:)./hrf_baseline;
        wgr_write_file(fname,dat,flag_nii_gii,v0)
        
        % write number of events
        fname = fullfile(outdir,[job.prefix,name,'_event_number',ext_nii_gii]);
        out.event_number{1} = fname;
        dat(smask_ind)=event_number;
        wgr_write_file(fname,dat,flag_nii_gii,v0)
    end
    if job.savedata.hrfmat_save
        save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'smask_ind','v0','-append')
    end
    
    
    if job.rmoutlier && flag_nii_gii %only for 3D data
        PARA_rm = PARA;
        pvalue_rm=job.para_global.pvalue_rm;
        Inpainted_method =job.para_global.Inpainted ;
        for i=1:3
            fname = fullfile(outdir,[job.prefix,name,'_Olrm_',HRF_para_str{i},ext_nii_gii]);
            eval(['out.Olrm_',HRF_para_str{i},'{1}=','''',fname,''';'])
            v0.fname = fname;
            dat(smask_ind)=PARA(i,:);
            if i==1
                id = (PARA_rm(1,:)<=0); %Be Careful
                PARA_rm(1,id)=nan;
                [PARA_rm(1,:),idx,outliers] = deleteoutliers(PARA_rm(1,:),pvalue_rm,1);
                fprintf('Response height outlier: negative/zeros #%d + extreme #%d, [%3.3f %3.3f];\n',nnz(id),length(idx),min(outliers),max(outliers));
                dat(smask_ind)=PARA_rm(1,:);
                dat=inpaint_nans3(dat,Inpainted_method);
                PARA_rm(1,:) = dat(smask_ind);
                id_rm = union(find(id),idx);
                dat2 = zeros(v0.dim);
                PARA_rm0=zeros(1,size(PARA,2));
                PARA_rm0(1,id_rm)=1;
                dat2(smask_ind)=PARA_rm0;
                v00 = v0;
                fname = fullfile(outdir,[job.prefix,name,'_outlier_NAN',ext_nii_gii]);
                out.outlier{1} = fname; 
                v00.fname = fname; 
                spm_write_vol(v00,dat2);  clear dat2
            else
                PARA_rm(i,id_rm)=nan;
                dat(smask_ind)=PARA_rm(i,:);
                dat=inpaint_nans3(dat,Inpainted_method);
                PARA_rm(i,:) = dat(smask_ind);
            end
            if job.savedata.hrfnii_save
                spm_write_vol(v0,dat);
            end
        end
        if job.savedata.hrfnii_save || job.savedata.hrfmat_save
            % write number of events
            fname = fullfile(outdir,[job.prefix,name,'_Olrm_event_number',ext_nii_gii]);
            out.Olrm_event_number{1} = fname;
            v0.fname = fname;
            event_number_rm  = event_number;
            event_number_rm(id_rm)=nan;
            dat(smask_ind)=event_number_rm;
            dat=inpaint_nans3(dat,Inpainted_method);
            event_number_rm = dat(smask_ind);
            if job.savedata.hrfnii_save
                spm_write_vol(v0,dat);
            end
            
            if flag_save_psc
                %write height PSC
                fname = fullfile(outdir,[job.prefix,name,'_Olrm_Height_PSC',ext_nii_gii]);
                out.Olrm_Height_PSC{1} = fname;
                v0.fname = fname;
                HeightPSC_rm  = PARA_rm(1,:)./hrf_baseline;
                HeightPSC_rm(HeightPSC_rm<=0)=nan;
                [HeightPSC_rm,idx,outliers] = deleteoutliers(HeightPSC_rm,pvalue_rm,1);
                fprintf('Response height (PSC) outlier: negative/zeros #%d + extreme #%d, [%3.3f %3.3f];\n',nnz(id),length(idx),min(outliers),max(outliers));
                dat(smask_ind)=HeightPSC_rm;
                dat=inpaint_nans3(dat,Inpainted_method);
                HeightPSC_rm = dat(smask_ind);
                if job.savedata.hrfnii_save
                    spm_write_vol(v0,dat);
                end
            end
            
            if job.savedata.hrfmat_save
                if flag_save_psc
                    save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'PARA_rm','event_number_rm','HeightPSC_rm','id_rm','-append')
                else
                    save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'PARA_rm','event_number_rm','id_rm','-append')
                end
            end
        end
    end
    
    if ~isempty(connroinfo) && any(conndata_flag~=2)
        conid = [conndata_flag~=2];
        data_nD = nan(Nscans, prod(v0.dim)); 
        data_nD(:,smask_ind) = data;
        wgr_conn_run(data_nD, connroinfo(conid,:),v0,name,outdir,flag_pval_pwgc,flag_nii_gii);    
        clear data_nD data
    end
 
 
    if job.savedata.deconv_save || any(conndata_flag~=1)
        % writing back deconvolved data into nifti/gifti file
        fname = fullfile(outdir,[job.prefix,name,ext_nii_gii]);
        out.deconv_data{1} = fname;
        if flag_nii_gii==1
            v1 = v; dat3 = zeros(v1(1).dim);
            for i=1:Nscans
                v1(i).fname = fname;
                v1(i).dt = [16,0]; 
                dat3(smask_ind) = data_deconv(i,:);
                if job.savedata.deconv_save
                    spm_write_vol(v1(i),dat3);
                end
            end
        else
            dat_surf= gifti;
            dat = zeros(nvar_data,Nscans);
            dat(smask_ind,:) = data_deconv';
            dat_surf.cdata = dat;
            save(dat_surf,fname,'Base64Binary');
        end
            
 
        if ~isempty(connroinfo) && any(conndata_flag~=1)
            conid =  [conndata_flag~=1];
            data_deconv3D = nan(Nscans, prod(v0.dim)); 
            data_deconv3D(:,smask_ind) = data_deconv;
            name2 = [name,'_deconv'];
            wgr_conn_run(data_deconv3D, connroinfo(conid,:),v0,name2,outdir,flag_pval_pwgc,flag_nii_gii);                            
        end
 
        if job.rmoutlier && job.para_global.rmoutlier_deconv
            v1 = v; dat3 = zeros(v1(1).dim);
            fname = fullfile(outdir,[job.prefix,name,'_Olrm',ext_nii_gii]);
            out.Olrm_deconv_data{1} = fname;
            for i=1:Nscans
                v1(i).fname = fname;
                v1(i).dt = [16,0]; 
                data_deconv_rm= data_deconv(i,:);
                data_deconv_rm(id_rm)=nan;
                dat3(smask_ind) = data_deconv_rm;
                dat3=inpaint_nans3(dat3,Inpainted_method);
                data_deconv3D(i,:) = dat3(:);
                if job.savedata.deconv_save
                    spm_write_vol(v1(i),dat3);
                end
            end
 
            if ~isempty(connroinfo)&& any(conndata_flag~=1)
                conid =  [conndata_flag~=1];
                name2 = [name,'_deconv_Olrm'];
                wgr_conn_run(data_deconv3D, connroinfo(conid,:),v0,name2,outdir,flag_pval_pwgc,flag_nii_gii);     
                clear data_deconv3D
            end
 
            clear data_deconv_rm
        end
        clear data_deconv data
    end
    clear dat
else %ROI-wise
 
    if ~isempty(job.connectivity)
        fprintf('Connectivity analysis...\n');
        [connroinfo,conndata_flag] = wgr_conn_check(job);
        if any(conndata_flag~=2)
            conid = [conndata_flag~=2];
            wgr_conn_run(data, connroinfo(conid,:),[],name,outdir,flag_pval_pwgc);  
        end
        if any(conndata_flag~=1)
            conid = [conndata_flag~=1];
            name2 = [name,'_deconv'];
            wgr_conn_run(data_deconv, connroinfo(conid,:),[],name2,outdir,flag_pval_pwgc);   
        end
    end
    
    if job.savedata.deconv_save
        try
            save(fullfile(outdir,[job.prefix,name, '.mat']),'data','data_deconv','-append');
        catch
            save(fullfile(outdir,[job.prefix,name, '.mat']),'data','data_deconv','-v7.3');
        end
    end
    clear data data_deconv
    out=[];
end
if job.savedata.job_save
    save(fullfile(outdir,[job.prefix,name, '_job.mat']), 'job');  
end

function wgr_ROI_sig_job(job)
ROI = job.Datasig; %cell file
[data_txt,mat_name]= wgr_check_ROIsig(ROI);
tmp = data_txt(:,1);
if job.para_global.combine_ROI
    data = cell2mat(tmp'); % combine all ROI together 
    tmp={}; tmp{1} = data;
    [~,outname,~] = fileparts(data_txt{1,2});
    fprintf('Combine all input ROIs for connectivity analysis, output is "*comb_%s*.mat"\n',outname)
end
nROI = length(tmp);
for i=1:nROI
    data = tmp{i};
    if isempty(data_txt{1,2})
        error('No data input!')
    else
        [fpath,name,~] = fileparts(data_txt{1,2});
        if job.para_global.combine_ROI
            job.raw_outname = ['combROI_',name];
        else
            if isempty(mat_name{i})
                job.raw_outname = name;
            else
                job.raw_outname = [name,'_',mat_name{i}];
            end
        end
    end
    outdir = job.outdir{1};
    if isempty(outdir)
        outdir = fpath;
    else
        if ~exist(outdir,'dir')
            mkdir(outdir)
        end
    end
    data = wgr_Denoising_job(job,data);
    flag_ROI = 1;
    
    if isfield(job,'HRFE') % deconvolution
        wgr_deconv_job(job,data,flag_ROI, outdir);
    else
        wgr_conn_job(job,data,flag_ROI, outdir);
    end
end

function [connroinfo,conndata_flag] = wgr_conn_check(job)
numBC = length(job.connectivity);
connroinfo = cell(numBC,5);
conndata_flag = zeros(numBC,1);
for j=1:numBC
    if isfield(job.connectivity{j},'conn')
        conns = job.connectivity{j}.conn;
    elseif isfield(job.connectivity{j},'connG')
        conns = job.connectivity{j}.connG;
    elseif isfield(job.connectivity{j},'connP')
        conns = job.connectivity{j}.connP;
    end
    if isfield(conns,'data4conn')
        conndata_flag(j,1) =conns.data4conn;
    end
    if isfield(conns,'genericROI')
        genericROI =conns.genericROI;
        if ~isfield(job,'ref_nii');
            job.ref_nii = job.images{1};
        end
        [mat,atlas,nii,atlasmeh,gii]= wgr_check_ROI(genericROI, job.ref_nii,job.para_global.delete_files);  
        connroinfo{j,1} = [mat;atlas;nii;atlasmeh;gii];
        connroinfo{j,2} =conns.Seed_ROI;
    else
        connroinfo{j,1} = [];
        connroinfo{j,2} = 1;
    end
    connroinfo{j,3} =conns.method;
    if isfield(conns,'morder')
        connroinfo{j,4} =conns.morder;
    end
    if isfield(conns,'ndinfo')
        connroinfo{j,5} =conns.ndinfo;
    end
    connroinfo{j,6} =conns.prefix;
end

function wgr_conn_job(job,data,flag_ROI, outdir, v, smask_ind, flag_nii_gii)
Nscans = size(data,1);
name = job.raw_outname;
flag_pval_pwgc = job.para_global.flag_pval_pwgc;
if ~flag_ROI % 2D/3D data
    if flag_nii_gii==1
        v0 = v(1);   v0.dt = [16,0]; v0.n = [1 1]; 
    else
        v0=v;
    end
    if ~isempty(job.connectivity)
        fprintf('Connectivity analysis...\n');
        connroinfo = wgr_conn_check(job);
    end
    if ~isempty(connroinfo) 
        data_nD = nan(Nscans, prod(v0.dim)); 
        data_nD(:,smask_ind) = data;
        wgr_conn_run(data_nD, connroinfo,v0,name,outdir,flag_pval_pwgc,flag_nii_gii);    
        clear data_nD data
    else
        warning('Please Configure Parameters for Connectivity Analysis')
    end
else %ROI-wise
    if ~isempty(job.connectivity)
        fprintf('Connectivity analysis...\n');
        connroinfo = wgr_conn_check(job);
        wgr_conn_run(data, connroinfo,[],name,outdir,1);
        clear data
    else
        warning('Please Configure Parameters for Connectivity Analysis')
    end
end

if job.job_save
    save(fullfile(outdir,[name, '_conn_job.mat']), 'job');  
end

function wgr_conn_run(data, connroinfo,v0,name,outdir,flag_pval_pwgc,flag_nii_gii);
%data: nobs x nvar (3D index)
fprintf('Connectivity analysis...\n ')
meastr = {'pwGC','CGC','Pearson','Spearman','PCGC'};
para_global = wgr_rsHRF_global_para;
regmode = para_global.regmode; % for GC
for j=1:size(connroinfo,1)
    roiid = connroinfo{j,1};  nroi = size(roiid,1);
    flag_seedROI  = connroinfo{j,2};
    conn_type = connroinfo{j,3};
    order = connroinfo{j,4};
    ndinfo = connroinfo{j,5}; 
    ROI = [];
    for i=1:nroi
        tmp = nanmean(data(:,roiid{i,1}),2) ;
        ROI(:,i) = tmp;
    end
    if nroi==0 %&& isempty(v0)
        ROI = data;
    end
    if ~flag_seedROI % seed map
        if flag_nii_gii==1
            ext_nii_gii = '.nii';
        else
            ext_nii_gii = '.gii';
        end
            
        smask_ind = find(var(data)>0);
        dat.data = data(:,smask_ind);
        dat.ROI = ROI;
        dat3 = zeros(v0.dim);
        if conn_type==1 || conn_type==2  || conn_type==5  
            if conn_type==2
                conn_type = 1; % pairwise       
                warning('Change Contional GC to Pairwise GC')
            elseif conn_type==5
                conn_type = 1; % pairwise       
                warning('Change Partially Contioned GC to Pairwise GC')
            end
            disp('Pairwise GC for seed based connectivity analysis')            
            [M] = wgr_GC(dat,order,ndinfo,conn_type,regmode,flag_pval_pwgc);
            ordeinfo = ['_order',num2str(order)];
            for i=1:nroi
                if nroi>1
                    tmp = [num2str(i),'_'];
                else
                    tmp = '';
                end
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_pwGC',ordeinfo,ext_nii_gii]);
                out.outflow_pwGC{i} = fname;
                gc = M.GC_Matrix_Out(i,:);
                dat3(smask_ind) = gc;
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
                
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_N_pwGC',ordeinfo,ext_nii_gii]);
                out.outflow_N_pwGC{i} = fname;
                dat3(smask_ind) = wgr_pwgc_2normal(gc,M.nobs,order);
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
                
                if flag_pval_pwgc
                    fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_pval_pwGC',ordeinfo,ext_nii_gii]);
                    out.outflow_pval_pwGC{i} = fname;
                    dat3(smask_ind) = M.pval_GC_Matrix_Out(i,:);
                    wgr_write_file(fname,dat3,flag_nii_gii,v0)
                end
                
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_pwGC',ordeinfo,ext_nii_gii]);
                out.inflow_pwGC{i} = fname;
                gc = M.GC_Matrix_In(i,:);
                dat3(smask_ind) = gc;
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
                
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_N_pwGC',ordeinfo,ext_nii_gii]);
                out.inflow_N_pwGC{i} = fname;
                dat3(smask_ind) = wgr_pwgc_2normal(gc,M.nobs,order);
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
                if flag_pval_pwgc
                    fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_pval_pwGC',ordeinfo,ext_nii_gii]);
                    out.inflow_pval_pwGC{i} = fname;
                    dat3(smask_ind) = M.pval_GC_Matrix_In(i,:);
                    wgr_write_file(fname,dat3,flag_nii_gii,v0)
                end
            end
        else
            [M] = wgr_FC(dat,conn_type);
            for i=1:nroi
                if nroi>1
                    tmp = [num2str(i),'_'];
                else
                    tmp = '';
                end
                
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_corr_',meastr{connroinfo{j,3}},ext_nii_gii]);
                out.corr{i} = fname;
                dat3(smask_ind) = M.Matrix_r(i,:);
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
                
                fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_Z_',meastr{connroinfo{j,3}},ext_nii_gii]);
                out.Z{i} = fname;
                dat3(smask_ind) = M.Matrix_z(i,:);
                wgr_write_file(fname,dat3,flag_nii_gii,v0)
            end
        end
        
    else %ROI to ROI
        dat.data = ROI;
        dat.ROI = [];
        if conn_type==1 || conn_type==2 || conn_type==5 % 1:pairwise  2:conditional 5: partically conditioned
            M = wgr_GC(dat,order,ndinfo, conn_type,regmode, flag_pval_pwgc);
            save(fullfile(outdir,[connroinfo{j,6},name,'_',meastr{connroinfo{j,3}},'.mat']),'M');
        else
            M = wgr_FC(dat,conn_type);
            save(fullfile(outdir,[connroinfo{j,6},name,'_Corr_',meastr{connroinfo{j,3}},'.mat']),'M');
        end
        
    end 
end

function wgr_write_file(fname,dat3,flag_nii_gii,v0)
if flag_nii_gii==1
    v0.fname = fname;
    spm_write_vol(v0,dat3);
else
    dat_surf= gifti;
    dat_surf.cdata = dat3;
    save(dat_surf,fname,'Base64Binary');
end

function c = wgr_pwgc_2normal(gc,nobs,order)
c = (nobs-order).*gc - (order-1)/3;
c(c<0)=0;
c = sqrt(c);

function [M] = wgr_GC(dat,order,ndinfo,flag_pw_cgc,regmode,flag_pval_pwgc,m);
% flag_pw_cgc, 1: pairwise GC,  2: conditional GC,  5: partially conditioned GC.
gcstr = {'Pairwise GC','Conditional GC','Pearson','Spearman','Partially Conditioned GC'};
data = dat.data;
[nobs,nvar] = size(data);
ROI =  dat.ROI;
[nobsROI,nROI] = size(ROI);
M.seed_num = nROI;
if nROI
    if nobsROI~=nobs
        error('different observations (ROI vs Data)')
    end
end
if nobs<10
    warning('Too few observations !')
end

fprintf('Data dimension: %d x %d\n',nobs,nvar);
if flag_pw_cgc==2 % CGC
    if nobs<nvar
        fprintf('#observation < #variable\n');
        error('CGC stop!')
    end
end

if nargin<7
    m = 5000; %block size
end

M.nvar=nvar;
M.nobs=nobs;
M.order = order;
M.GC_type = gcstr{flag_pw_cgc};
M.ndinfo = ndinfo;
if flag_pw_cgc==5 % if nnz(ndinfo)==1
    nd = ndinfo(1);
    if isnan(nd)
        nd = 6;
    end
    ndmax = ndinfo(2); 
    if isnan(ndmax)
        ndmax = nd+1;
    end
end

nbin = ceil(nvar/m);
indX={}; indY={};
matrix={};
matrix_out = {};
matrix_in={};
p_matrix_out={};
p_matrix_in={};

if nROI %ROI to data
    nbinROI = ceil(nROI/m);
    for i=1:nbinROI
        if i~=nbinROI
            ind_X = (i-1)*m+1:i*m ;   
        else
            ind_X = (i-1)*m+1:nROI ;
        end
        indX{i} = ind_X; 
        for j=1:nbin
            if j~=nbin
                ind_Y = (j-1)*m+1:j*m ;
            else
                ind_Y = (j-1)*m+1:nvar ;
            end
            indY{j} = ind_Y; 
            if flag_pval_pwgc
                [matrix_out{i,j}, matrix_in{i,j}, p_matrix_out{i,j}, p_matrix_in{i,j}] = ...
                    wgr_seedGC(ind_X,ind_Y,ROI,data,order,1,regmode); 
            else
                [matrix_out{i,j}, matrix_in{i,j}] = wgr_seedGC(ind_X,ind_Y,ROI,data,order,0,regmode); %only pairwise.
            end
        end
    end
else % data to data
    if flag_pw_cgc==5 
        [y_inform_gain, cind] = wgr_init_partial_conditioning(data,[],ndmax,order);
        M.information_gain = y_inform_gain;
        M.condition_id = cind;
    end
    if flag_pw_cgc==1
        [M.GC_Matrix, M.pval_Matrix] = wgr_PWGC(data,order,regmode,flag_pval_pwgc);
    end
    if flag_pw_cgc==2
        [M.GC_Matrix, M.pval_Matrix] = rsHRF_mvgc(data',order,regmode,0,flag_pval_pwgc);
    end
    
    if flag_pw_cgc==5
        [M.GC_Matrix, M.pval_Matrix] = wgr_PCGC(data,order,cind,nd,regmode,flag_pval_pwgc);
    end
end

if nROI
    M.GC_Matrix_Out = nan(nROI,nvar);
    M.GC_Matrix_In = M.GC_Matrix_Out ;
    for i=1:nbinROI
        for j=1:nbin
            M.GC_Matrix_Out(indX{i},indY{j}) =  matrix_out{i,j};
            M.GC_Matrix_In(indX{i},indY{j}) =  matrix_in{i,j};
            if flag_pval_pwgc
                M.pval_GC_Matrix_Out(indX{i},indY{j}) =  p_matrix_out{i,j};
                M.pval_GC_Matrix_In(indX{i},indY{j}) =  p_matrix_in{i,j};
            end
        end
    end 

else
    
    if flag_pw_cgc==1
        M.GC_Matrix_N = wgr_pwgc_2normal(M.GC_Matrix,nobs,order);
    end        
end
clear matrix* p_* indX indY

function [F,pvalue] = wgr_PWGC(data,order,regmode,flag_pval);

[nvar] = size(data,2);
F = zeros(nvar);
if flag_pval
    pvalue = nan(nvar);
end
for drive=1:nvar
    for target=1:nvar
        if drive~=target
            dat = data(:,[drive target])';
            [F0,p0] = rsHRF_mvgc(dat,order,regmode,0,flag_pval);
            F(drive,target) = F0(1,2);
            pvalue(drive,target)  = p0(1,2);
            F(target,drive) = F0(2,1);
            pvalue(target,drive)  = p0(2,1);
        end
    end
end

function [F,pvalue] = wgr_PCGC(data,order,cind,nd,regmode,flag_pval);
[nvar] = size(data,2);
F = zeros(nvar);
if flag_pval
    pvalue = nan(nvar);
end
for drive=1:nvar
    for target=1:nvar
        if drive~=target
            zid = setdiff(cind(drive,:),target,'stable');
            dat = data(:,[drive target zid(1:nd)])';
            [F(drive,target),pvalue(drive,target)] = ...
                rsHRF_mvgc(dat,order,regmode,1,flag_pval);
        end
    end
end

function [gc_out, gc_in, p_out, p_in] = wgr_seedGC(ind_X,ind_Y,ROI,data1,order,flag_pval,regmode);

[nvar1] = length(ind_X);
[nvar2] = length(ind_Y);
gc_out = zeros(nvar1,nvar2);
gc_in = gc_out;
if flag_pval
    p_out = nan(nvar1,nvar2);
    p_in = p_out;
end
for drive=1:nvar1
    for target=1:nvar2
        data = [ROI(:,ind_X(drive)) data1(:,ind_Y(target))]';
        [F,pvalue] = rsHRF_mvgc(data,order,regmode,0,flag_pval);
        if length(F)>1
            gc_out(drive,target) = F(1,2); 
            gc_in(drive,target)  = F(2,1); 

            p_out(drive,target)  = pvalue(1,2);
            p_in(drive,target)   = pvalue(2,1); 
        end
    end
end

function [M] = wgr_FC(dat,conn_type,m);
% conn_type, 3: pearson,  4: conn_type

if conn_type==3
    con_type =  'Pearson';
elseif conn_type==4
    con_type = 'Spearman';
end

data = dat.data;
[nobs,nvar] = size(data);
ROI =  dat.ROI;
[nobsROI,nROI] = size(ROI);
M.seed_num = nROI;
if nROI
    if nobsROI~=nobs
        error('different observations (ROI vs Data)')
    end
end
if nobs<10
    warning('Too few observations !')
end

fprintf('Data dimension: %d x %d\n',nobs,nvar);

if nargin<3
    m = 5000; %block size
end

M.nvar=nvar;
M.nobs=nobs;

nbin = ceil(nvar/m);
indX={}; indY={};
matrix={};

if nROI %ROI to data
    nbinROI = ceil(nROI/m);
    for i=1:nbinROI
        if i~=nbinROI
            ind_X = (i-1)*m+1:i*m ;   
        else
            ind_X = (i-1)*m+1:nROI ;
        end
        indX{i} = ind_X; 
        for j=1:nbin
            if j~=nbin
                ind_Y = (j-1)*m+1:j*m ;
            else
                ind_Y = (j-1)*m+1:nvar ;
            end
            indY{j} = ind_Y;  
            matrix_r{i,j} = corr(ROI(:,ind_X),data(:,ind_Y), 'type', con_type); 
            matrix_z{i,j} = atanh(matrix_r{i,j}) ;
        end
    end
else % data to data
    for i=1:nbin
        if i~=nbin
            ind_X = (i-1)*m+1:i*m ;   
        else
            ind_X = (i-1)*m+1:nvar ;
        end
        indX{i} = ind_X; 
        for j=1:nbin
            if j~=nbin
                ind_Y = (j-1)*m+1:j*m ;
            else
                ind_Y = (j-1)*m+1:nvar ;
            end
            indY{j} = ind_Y; 
            [matrix_r{i,j},matrix_p{i,j}] = corr(data(:,ind_X),data(:,ind_Y), 'type',con_type); 
            matrix_z{i,j} = atanh(matrix_r{i,j}) ;
        end
    end
end

if nROI
    M.Matrix_r = nan(nROI,nvar);
    M.Matrix_z = M.Matrix_r;
    for i=1:nbinROI
        for j=1:nbin
            M.Matrix_r(indX{i},indY{j}) =  matrix_r{i,j};
            M.Matrix_z(indX{i},indY{j}) =  matrix_z{i,j};
        end
    end 
    
else
    M.Matrix_r = nan(nvar,nvar);
    for i=1:nbin
        for j=1:nbin
            M.Matrix_r(indX{i},indY{j}) =  matrix_r{i,j};
            M.Matrix_z(indX{i},indY{j}) =  matrix_z{i,j};
            M.Matrix_pval(indX{i},indY{j}) =  matrix_p{i,j};
        end
    end 
    M.Matrix_r(1:nvar+1:end)=0; %remove diag value.
    M.Matrix_z(1:nvar+1:end)=0; %remove diag value.
    M.Matrix_pval(1:nvar+1:end)=1; %remove diag value.
end
clear matrix* indX indY

function [y, ind] = wgr_init_partial_conditioning(data,seed_signal,ndmax,order)
[N,nvar] = size(data);
X=cell(nvar,1);
past_ind = repmat([1:order],N-order,1) + repmat([0:N-order-1]',1,order);
for i=1:nvar
    past_data= reshape(data(past_ind,i),N-order,order);
    X{i}= past_data - repmat(mean(past_data), N-order, 1); %%remove mean
end

ind=zeros(nvar,ndmax);
y=ind;

if ~isempty(seed_signal)
   [N2,nvar2] = size(seed_signal);
   if N~=N2
       error('check #observation of seed_signal !')
   else
       ind = nan(nvar2,ndmax);  y = ind;
       parfor drive=1:nvar2
           past_data = reshape(seed_signal(past_ind,drive),N-order,order);
           drive_sig = past_data - repmat(mean(past_data), N-order, 1);
           [y(drive,:), ind(drive,:)]=wgr_information_gain(drive_sig,X,nvar,ndmax);
       end
   end
else
    ind = nan(nvar,ndmax);  y = ind;
    parfor drive=1:nvar
        [y(drive,:), ind(drive,:)]=wgr_information_gain(drive,X,nvar,ndmax); %do not allow to dynamic plot
    end
end

function [y, ind]= wgr_information_gain(drive,X,nvar,ndmax)
if nnz(drive)==1 
    indici=setdiff(1:nvar,drive);%,'stable');
    t=X{drive};
else %seed signal
    indici= 1:nvar ;
    t=drive;
end

Zt=[]; ind=nan(1,ndmax); y = ind;
for nd=1:ndmax
    n1=length(indici);
    z=zeros(n1,1);
    for k=1:n1
        Zd=[Zt X{indici(k)}];
        z(k)= wgr_MI_gaussian(Zd,t);
    end
    [y(1,nd), id] = max(z);
    Zt = [Zt, X{indici(id)}];
    ind(1,nd) = indici(id);
    indici = setdiff(indici, indici(id));%,'stable');
end
 
function mi = wgr_MI_gaussian(Y,X)
% X: N x k_x  Y: N x k_y 
% condtional mutual information
% mi = gaussian mutual information X,Y
XY_cov = X'*Y; %%here we delete the factor: 1/N.
X_cov = X'*X;
Y_cov = Y'*Y;
all_cov = [X_cov XY_cov;XY_cov' Y_cov];
mi = 0.5*log(abs(det(X_cov)*det(Y_cov)/det(all_cov)));
return
