function rsHRF(job,flag)
%% resting state BOLD-fMRI HRF deconvolution
if nargin>0
    para_global = wgr_rsHRF_global_para;
    switch flag;
        case 'vox'
            flag_delete = para_global.delete_files; % delete temporary files (generated wm/csf/brainmask) 
            TR = job.HRFE.TR ; 
            nii_file = job.images; %cell file
            nii_mask = job.mask{1};
            if isempty(nii_mask)
                nii_mask = para_global.mask_nii ;
            end
            
            data = [];
            fprintf('Reading NIFTI Data ...\n')
            Nscans = length(nii_file);
            for i=Nscans:-1:1
                v(i) = spm_vol(nii_file{i});
                tmp  = spm_read_vols(v(i));
                data(i,:) = tmp(:);
            end       
            
            outdir = job.outdir{1};
            if isempty(outdir)
                outdir = fileparts(v(1).fname);
            end
            
            if  isfield(job,'genericROI')    
                ROI = job.genericROI;
            else
                ROI = [];
            end
            
            if ~isempty(ROI);
                flag_ROI = 1;
                [mat,atlas,nii]= wgr_check_ROI(ROI, nii_file{1},flag_delete);                
                roinfo = [mat;atlas;nii];
                nvar = size(roinfo,1);
                data_ROI = zeros(Nscans,nvar);
                for j=1:nvar
                    data_ROI(:,j) = mean(data(:,roinfo{j,1}),2);
                end
            else
                flag_ROI=0;
            end
            
            if ~isempty(nii_mask);
                nii_masknew  = spm_file(nii_mask,'prefix','bin_');
                spm_imcalc({nii_file{1};nii_mask},nii_masknew,'(i1~=0).*(i2~=0)',{0,0,0,16});
                v0 = spm_vol(nii_masknew);
                datam = spm_read_vols(v0);
                voxel_ind = find(datam); clear datam
                if flag_delete
                    delete(v0.fname)
                end
            else
                data_var = nanvar(data,0,1);
                voxel_ind = find(data_var); clear data_var
            end
            data_gm = mean(data(:,voxel_ind),2);
            fprintf('# %d voxels inside defaut mask ...\n',length(voxel_ind))
            
            covariates = job.Denoising.generic;
            fprintf('Reading Covariates ...\n')
            if ~isempty(covariates)
                [txt,mat,nii]= wgr_check_covariates(covariates);     
                if ~isempty(nii)
                    nii2 = {}; 
                    for i=1:length(nii)
                        nii2{1,i} = spm_file(nii{i,1},'prefix','bin_');
                        spm_imcalc({nii_file{1};nii{i,1}},nii2{1,i},'(i1~=0).*(i2>0.9)',{0,0,0,16});
                        v0 = spm_vol(nii2{1,i});
                        datam = spm_read_vols(v0);
                        if flag_delete
                            delete(v0.fname)
                        end
                        id = find(datam==1); clear datam
                        if isempty(id)
                            error(['no survived voxels (>0.9) in ',nii{i,1}])
                        else
                            if para_global.aCompcor,                                 
                                numcomps = para_global.aCompcor_numcomps;
                                z = data(:,id);
                                z(isnan(z))=0;
                                m = size(z,2) ;
                                if m<numcomps
                                    error(['# components = ',num2str(m),', < ', num2str(numcomps)])
                                else
                                    [pc,score,latent] = princomp(z,'econ');
                                end
                                nii2{1,i} = z*pc(:,1:numcomps);
                            else
                                nii2{1,i} = mean(data(:,id),2);
                            end
                        end
                    end
                    nii = cell2mat(nii2);
                    clear nii2 datam
                else
                    nii = [];
                end
                if para_global.aCompcor
                    nuisance = [txt mat nii];
                else ~para_global.aCompcor & para_global.globa_reg
                    nuisance = [txt mat nii data_gm];
                end
            else
                nuisance = [];
            end
            %% 
            if flag_ROI
                data = data_ROI;
            else
                data = data(:,voxel_ind);
            end
            if length(job.Denoising.bands)>1
                error('please check your setting for band-pass filtering \n')
            elseif length(job.Denoising.bands)==1
                Bands = job.Denoising.bands{1};   
            else
                Bands =[];
            end
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
            
            fprintf('Removing regressors & Detrend................. \n ')
            nuisance = [nuisance_detrend nuisance];
            if ~isempty(nuisance)
                nuisance = bsxfun(@minus, nuisance, sum(nuisance, 1)./Nscans);% Remove the mean
                nuisance = [orth(nuisance) ones(Nscans, 1)];
                mu   = sum(data, 1)./Nscans;    % mean 
                data = data - nuisance*(nuisance\data);
                data = bsxfun(@plus, data, mu);
            end
            
            if ~isempty(Bands)
                fprintf('Frequency filter.................... \n ')
                data = wgr_band_filter(data, TR,Bands);
            end
            
            if job.Denoising.Despiking;
                fprintf(' Temporal despiking with a hyperbolic tangent squashing function... \n')
                data = wgr_despiking(data);
            end           
            
            %% HRF deconvolution
            para.TR = TR;
            %%% the following parameter (upsample grid) can be > 1 only for Canonical.
            %%% Set = 1 for FIR
            para.T  = job.HRFE.fmri_t; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid
            para.T0 = job.HRFE.fmri_t; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
            if para.T==1
                para.T0 = 1;
            end
            min_onset_search = min(job.HRFE.mdelay); % minimum delay allowed between event and HRF onset (seconds)
            max_onset_search = max(job.HRFE.mdelay); % maximum delay allowed between event and HRF onset (seconds)
            para.dt  = para.TR/para.T; % fine scale time resolution.
            if job.HRFE.hrfm==1
                para.TD_DD = 2; % time and dispersion derivative
            elseif job.HRFE.hrfm==3
                para.TD_DD = 1; % time and dispersion derivative
            end
            para.AR_lag = job.HRFE.cvi; % AR(1) noise autocorrelation.
            para.thr = job.HRFE.thr; % (mean+) para.thr*standard deviation threshold to detect event.
            para.len = job.HRFE.hrflen; % length of HRF, in seconds
            para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
            para.localK = para_global.localK;
            temporal_mask = logical(ones(Nscans,1));
            if job.HRFE.hrfm==1 | job.HRFE.hrfm==3
                tic
                [beta_hrf, bf, event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data,para, temporal_mask);
                hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
            elseif job.HRFE.hrfm==2 | job.HRFE.hrfm==4
                tic
                if  job.HRFE.hrfm==2
                    para.estimation = 'FIR';
                elseif job.HRFE.hrfm==4
                    para.estimation = 'sFIR';
                end
                para.T=1; % this needs to be = 1 for FIR
                [hrfa,event_bold] = wgr_rsHRF_FIR(data,para, temporal_mask);
            end
            
            nvar = size(hrfa,2);
            PARA = zeros(3,nvar);
            for voxel_id=1:nvar
                hrf1 = hrfa(:,voxel_id);
                [PARA(:,voxel_id)] = wgr_get_parameters(hrf1,para.TR/para.T);% estimate HRF parameter
            end
            fprintf('\nDone HRF estimation %8.2f seconds\n',toc)
            
            fprintf('Deconvolving HRF ...\n');
            tic
            T = round(para.len/TR);
            if para.T>1
                hrfa_TR = resample(hrfa,1,para.T);
            else
                hrfa_TR = hrfa;
            end
            for voxel_id=1:nvar
                hrf=hrfa_TR(:,voxel_id);
                H=fft([hrf; zeros(Nscans-length(hrf),1)]);
                M=fft(data(:,voxel_id));
                data_deconv(:,voxel_id) = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));
                event_number(voxel_id)=length(event_bold{1,voxel_id});
            end
            toc
            fprintf('Save HRF deconvolution results....\n');
            v0 = v(1);   v0.dt = [16,0]; v0.n = [1 1];
            [fpath,name,ext] = fileparts(v0.fname);
            save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'para', 'hrfa', 'event_bold', 'event_number','PARA','-v7.3');
            
            if ~flag_ROI
                %write nifti maps of the three HRF parameters
                HRF_para_str = {'Height.nii', 'Time2peak.nii','FWHM.nii'};
                dat= zeros(v0.dim);
                for i=1:3
                    v0.fname = fullfile(outdir,[job.prefix,name,'_',HRF_para_str{i}]);
                    dat(voxel_ind)=PARA(i,:);
                    spm_write_vol(v0,dat);
                end
                % write number of events
                v0.fname = fullfile(outdir,[job.prefix,name,'_event_number.nii']);
                dat(voxel_ind)=event_number;
                spm_write_vol(v0,dat);
                
                if job.rmoutlier
                    PARA_rm = PARA;
                    pvalue_rm=para_global.pvalue_rm;
                    Inpainted_method =para_global.Inpainted ;
                    for i=1:3
                        v0.fname = fullfile(outdir,[job.prefix,name,'_Olrm_',HRF_para_str{i}]);
                        dat(voxel_ind)=PARA(i,:);
                        if i==1
                            id = (PARA_rm(1,:)<=0);
                            PARA_rm(1,id)=nan;
                            [PARA_rm(1,:),idx,outliers] = deleteoutliers(PARA_rm(1,:),pvalue_rm,1);
                            fprintf('height outlier: negative/zeros #%d + extreme #%d, [%3.3f %3.3f];\n',nnz(id),length(idx),min(outliers),max(outliers));
                            dat(voxel_ind)=PARA_rm(1,:);
                            dat=inpaint_nans3(dat,Inpainted_method);
                            PARA_rm(1,:) = dat(voxel_ind);
                            id_rm = union(find(id),idx);
                            dat2 = zeros(v0.dim);
                            PARA_rm0=zeros(1,size(PARA,2));
                            PARA_rm0(1,id_rm)=1;
                            dat2(voxel_ind)=PARA_rm0;
                            v00 = v0;
                            v00.fname = fullfile(outdir,[job.prefix,name,'_outlier_NAN.nii']);
                            spm_write_vol(v00,dat2);  clear dat2
                        else
                            PARA_rm(i,id_rm)=nan;
                            dat(voxel_ind)=PARA_rm(i,:);
                            dat=inpaint_nans3(dat,Inpainted_method);
                            PARA_rm(i,:) = dat(voxel_ind);
                        end
                        spm_write_vol(v0,dat);
                    end
                    % write number of events
                    v0.fname = fullfile(outdir,[job.prefix,name,'_Olrm_event_number.nii']);
                    dat(voxel_ind)=event_number;
                    event_number_rm  = event_number;
                    event_number_rm(1,id_rm)=nan;
                    dat(voxel_ind)=event_number_rm;
                    dat=inpaint_nans3(dat,Inpainted_method);
                    event_number_rm = dat(voxel_ind);
                    spm_write_vol(v0,dat);
                    save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'PARA_rm','event_number_rm','id_rm','-append')
                end

                % writing back deconvolved data into nifti file
                v1 = v; dat3 = zeros(v1(1).dim);
                for i=1:Nscans
                    [fpath,name,ext] = fileparts(v1(i).fname);
                    v1(i).fname = fullfile(outdir,[job.prefix,name,ext]);
                    v1(i).dt = [16,0]; 
                    dat3(voxel_ind) = data_deconv(i,:);
                    spm_write_vol(v1(i),dat3);
                end
                if job.rmoutlier & para_global.rmoutlier_deconv
                    v1 = v; dat3 = zeros(v1(1).dim);
                    
                    for i=1:Nscans
                        [fpath,name,ext] = fileparts(v1(i).fname);
                        v1(i).fname = fullfile(outdir,[job.prefix,name,'_Olrm',ext]);
                        v1(i).dt = [16,0]; 
                        data_deconv_rm= data_deconv(i,:);
                        data_deconv_rm(id_rm)=nan;
                        dat3(voxel_ind) = data_deconv_rm;
                        dat3=inpaint_nans3(dat3,Inpainted_method);
                        spm_write_vol(v1(i),dat3);
                    end
                    clear data_deconv_rm
                end
                clear data_deconv
            else
                save(fullfile(outdir,[job.prefix,name, '.mat']), 'roinfo','data','data_deconv','-v7.3');
                clear data data_deconv
            end
            save(fullfile(outdir,[job.prefix,name, '_job.mat']), 'job','para_global');            
        case 'sig'
            wgr_ROI_sig_deconvolution(job,para_global)
    end
    return
else
    code_path = fileparts(which('rsHRF.m')) ;
    a = imread(fullfile(code_path,'rsHRF_logo.png'));
    S.fig = figure('Visible','on',...
        'numbertitle','off',...    
        'menubar','none',...         
        'units','normalized',...
        'color','w',...
        'position',[0.4   0.4    0.15    0.3],...[563    98   480   360],...
        'name',['rsHRF v1.0(',getenv('USERNAME'),')'],...
        'resize','off');
    
    axes('parent',S.fig,'units','normalized','position',[0 0.59 1 0.43]);
    imagesc(a);
    axis off; 
    text(80,360,'Resting State HRF','Color',[1 0.6 0],'Fontsize',13,'Fontweight','bold');

    
    %% ROI signal HRF deconvolution
    S.pb(1) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.1,0.05,0.8,0.15],...
        'style','pushbutton',...'pushbutton',...
        'string','Signals',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.HRF.sig_rsHRF'');',...
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.35,... 
        'fontweight','bold');


    %% ROI-wise HRF deconvolution
    S.pb(2) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.1,0.25,0.8,0.15],...
        'style','pushbutton',...'pushbutton',...
        'string','ROIs',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.HRF.ROI_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.35,... 
        'fontweight','bold');
    
    %% voxelwise HRF deconvolution
    S.pb(3) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.1,0.45,0.8,0.15],...
        'style','pushbutton',...'pushbutton',...
        'string','Voxels',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.HRF.vox_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.35,... 
        'fontweight','bold');
                                                                                                                                                      
end


function x1 = wgr_despiking(x1)
my=repmat(median(x1,1),[size(x1,1),1]);
sy=repmat(4*median(abs(x1-my)),[size(x1,1),1]);
x1=my+sy.*tanh((x1-my)./max(eps,sy));

function  x = wgr_band_filter(x, TR,Bands)
x = conn_filter(TR,Bands,x,'full') +...%mean-center by default
     repmat(mean(x,1),size(x,1),1);    %mean back in

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
        case 'full',
            idx=find(f<filter(1)*(rt*size(fy,1))|f>=filter(2)*(rt*size(fy,1)));
            %idx=idx(idx>1);
            fy(idx,:)=0;
            k=1; %2*size(fy,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(fy,1)-numel(idx));
            y=real(ifft(fy,[],1))*k;
        case 'partial',
            idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
            %if ~any(idx==1), idx=[1,idx]; end
            y=fy(idx,:);
    end
end


function [mat,atlas,nii]= wgr_check_ROI(ROI, funii, flag_delete)
mat = {};  atlas = {}; nii  = {};
k1 =1 ; k2 = 1; k3 = 1;
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
            mat{j,1} = id;
            mat{j,2} = tmp2;
            mat{j,3} = '';
            continue;
        end
    end
        
    flag = isfield(tmp,'atlas');
    if flag
        for j=1: length(tmp.atlas)
            tmp2 = strcat(tmp.atlas{j});
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
        continue;
    end
    
    flag = isfield(tmp,'images');
    if flag
        for j=1: length(tmp.images)
            tmp2 = strcat(tmp.images{j});            
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
        continue;
    end
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
            [fpath,name,ext] = fileparts(tmp2);
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
        
    flag = isfield(tmp,'images');
    if flag
        for j=1: length(tmp.images)
            tmp2 = strcat(tmp.images{j});
            nii{k3,1} = tmp2;
            k3 = k3+1;
        end
        continue;
    end
end
txt = cell2mat(txt);
mat = cell2mat(mat);

function wgr_ROI_sig_deconvolution(job,para_global)
TR = job.HRFE.TR ; 
ROI = job.Datasig; %cell file
[mat_name,mat_txt]= wgr_check_ROIsig(ROI);
tmp = mat_txt(:,1)';
data = cell2mat(tmp);
Nscans = size(data,1);
datainfo = mat_txt(:,2);

if isempty(mat_txt{1,2})
    error('No data input!')
else
    [fpath,name,ext] = fileparts(mat_txt{1,2});
end

outdir = job.outdir{1};
if isempty(outdir)
    outdir = fpath;
end
covariates = job.Denoising.generic;
fprintf('Reading Covariates ...\n')
if ~isempty(covariates)
    [txt,mat,nii]= wgr_check_covariates(covariates);     
    nuisance = [txt mat];
else
    nuisance = [];
end

if length(job.Denoising.bands)>1
    error('please check your setting for band-pass filtering \n')
elseif length(job.Denoising.bands)==1
    Bands = job.Denoising.bands{1};   
else
    Bands =[];
end
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
            
fprintf('Removing regressors & Detrend................. \n ')
nuisance = [nuisance_detrend nuisance];
if ~isempty(nuisance)
    nuisance = bsxfun(@minus, nuisance, sum(nuisance, 1)./Nscans);% Remove the mean
    nuisance = [orth(nuisance) ones(Nscans, 1)];
    mu   = sum(data, 1)./Nscans;    % mean 
    data = data - nuisance*(nuisance\data);
    data = bsxfun(@plus, data, mu);
end
            
if ~isempty(Bands)
    fprintf('Frequency filter.................... \n ')
    data = wgr_band_filter(data, TR,Bands);
end
            
if job.Denoising.Despiking;
    fprintf(' Temporal despiking with a hyperbolic tangent squashing function... \n')
    data = wgr_despiking(data);
end           
            
%% HRF deconvolution
para.TR = TR;
%%% the following parameter (upsample grid) can be > 1 only for Canonical.
%%% Set = 1 for FIR
para.T  = job.HRFE.fmri_t; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid
para.T0 = job.HRFE.fmri_t; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
if para.T==1
    para.T0 = 1;
end
min_onset_search = min(job.HRFE.mdelay); % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = max(job.HRFE.mdelay); % maximum delay allowed between event and HRF onset (seconds)
para.dt  = para.TR/para.T; % fine scale time resolution.
if job.HRFE.hrfm==1
    para.TD_DD = 2; % time and dispersion derivative
elseif job.HRFE.hrfm==3
    para.TD_DD = 1; % time and dispersion derivative
end
para.AR_lag = job.HRFE.cvi; % AR(1) noise autocorrelation.
para.thr = job.HRFE.thr; % (mean+) para.thr*standard deviation threshold to detect event.
para.len = job.HRFE.hrflen; % length of HRF, in seconds
para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
para.localK = para_global.localK;
temporal_mask = logical(ones(Nscans,1));
if job.HRFE.hrfm==1 | job.HRFE.hrfm==3
    tic
    [beta_hrf, bf, event_bold] = wgr_rsHRF_estimation_canonhrf2dd_par2(data,para, temporal_mask);
    hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
elseif job.HRFE.hrfm==2 | job.HRFE.hrfm==4
    tic
    if  job.HRFE.hrfm==2
        para.estimation = 'FIR';
    elseif job.HRFE.hrfm==4
        para.estimation = 'sFIR';
    end
    para.T=1; % this needs to be = 1 for FIR
    [hrfa,event_bold] = wgr_rsHRF_FIR(data,para, temporal_mask);
end

nvar = size(hrfa,2);
PARA = zeros(3,nvar);
for voxel_id=1:nvar
    hrf1 = hrfa(:,voxel_id);
    [PARA(:,voxel_id)] = wgr_get_parameters(hrf1,para.TR/para.T);% estimate HRF parameter
end
fprintf('Done HRF estimation %8.2f seconds\n',toc)

fprintf('Deconvolving HRF ...\n');
tic
T = round(para.len/TR);
if para.T>1
    hrfa_TR = resample(hrfa,1,para.T);
else
    hrfa_TR = hrfa;
end
for voxel_id=1:nvar
    hrf=hrfa_TR(:,voxel_id);
    H=fft([hrf; zeros(Nscans-length(hrf),1)]);
    M=fft(data(:,voxel_id));
    data_deconv(:,voxel_id) = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));
    event_number(voxel_id)=length(event_bold{1,voxel_id});
end
toc
fprintf('Save HRF deconvolution results....\n');
save(fullfile(outdir,[job.prefix,name, 'hrf.mat']), 'para', 'hrfa', 'event_bold', 'event_number','PARA','-v7.3');
save(fullfile(outdir,[job.prefix,name, '_deconv.mat']), 'datainfo','data','data_deconv','-v7.3');
save(fullfile(outdir,[job.prefix,name, 'job.mat']), 'job');

function [mat_name,mat_txt]= wgr_check_ROIsig(ROI)
mat_txt = {}; mat_name  = {};
k1 =1 ; k2 = 1; k3 = 1;
for i=1:length(ROI)
    tmp  = ROI{i};
    flag = isfield(tmp,'name');
    if flag
        if size(tmp.name,1)>1
            error('more than 1 variable names')
        else
            tmp2 = strcat(tmp.name);
            mat_name = tmp2;            
        end
        break
    end
end

for i=1:length(ROI)
    tmp  = ROI{i};
    flag = isfield(tmp,'sigdata');
    if flag
        for j=1: size(tmp.sigdata,1)
            tmp2 = tmp.sigdata{j};
            [fpath,name,ext] = fileparts(tmp2);
            if strfind(ext,'txt')
                mat_txt{j,1} = load(tmp2);
                mat_txt{j,2} = tmp2;
            elseif strfind(ext,'mat')
                tmp3 = load(tmp2);
                try
                    eval(['mat_txt{j,1} = tmp3.',mat_name,';']);
                catch
                    tmp4 = fieldnames(tmp3);
                    if length(tmp4)==1
                        eval(['mat_txt{j,1} = tmp3.',tmp4{1},';']);
                    else
                        error(['there are more than one variable in mat file, and no variable ''',mat_name,'" in mat file']);
                    end
                end
                mat_txt{j,2} = tmp2;
            end
            continue;
        end
    end  
end