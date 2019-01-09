function rsHRF(varargin)
%% resting state BOLD-fMRI HRF deconvolution and connectivity analysis
if nargin==2
    job = varargin{1};  flag = varargin{2};  
    if ~isfield(job,'para_global')
        para_global = wgr_rsHRF_global_para;
        job.para_global = para_global; clear para_global
    end
    switch flag
        case 'vox'
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
           
            if isfield(job,'HRFE') % deconvolution
                wgr_deconv_job(job,data,flag_ROI, outdir, v, voxel_ind)
            else
                wgr_conn_job(job,data,flag_ROI, outdir, v, voxel_ind)
            end
        case 'sig'
            wgr_ROI_sig_job(job)         
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
        'position',[0.4   0.3    0.15    0.38],...[563    98   480   360],...
        'name',['rsHRF v2.0(',getenv('USERNAME'),')'],...
        'resize','off');
    
    axes('parent',S.fig,'units','normalized','position',[0 0.59 1 0.43]);
    imagesc(a);
    axis off; 
    
    %% ROI signal HRF deconvolution
    S.pb(1) = uicontrol('parent',S.fig,...
        'unit','norm',...
        'pos',[0.1,0.05,0.8,0.15],...
        'style','pushbutton',...'pushbutton',...
        'string','Signals',...
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.sig_rsHRF'');',...
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
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.ROI_rsHRF'');',...        
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
        'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.vox_rsHRF'');',...        
        'backgroundc','w',...
        'foregroundcolor',0*[1 1 1],...     'fontname','Segoe Script',...
        'fontname','Calibri',...'Times New Roman',...
        'fontunits', 'normalized',... 
        'fontsize', 0.35,... 
        'fontweight','bold');
    
    if nargin==1
        Modality=upper(varargin{1});
        switch Modality
            case 'CONN'  
                set(S.pb(1),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.sig_conn'');'); 
                set(S.pb(2),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.ROI_conn'');'); 
                set(S.pb(3),'CallBack','spm_jobman(''interactive'','''',''spm.tools.rsHRF.vox_conn'');');
                set(S.pb(:),'foregroundcolor',[0.3 0.6 0])
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
                        nii2{1,i} = mean(data(:,id),2);
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
    data = wgr_band_filter(data, TR,Bands);
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

function wgr_deconv_job(job,data,flag_ROI, outdir, v, voxel_ind)
Nscans = size(data,1);
name = job.raw_outname;
flag_pval_pwgc = job.para_global.flag_pval_pwgc;

%% HRF deconvolution
para.TR = job.HRFE.TR ; 
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
para.localK = job.para_global.localK;
temporal_mask = true(Nscans,1); % temporal mask.
if job.HRFE.hrfm==1 || job.HRFE.hrfm==3
    tic
    [beta_hrf, bf, event_bold] = wgr_rshrf_estimation_canonhrf2dd_par2(data,para, temporal_mask);
    hrfa = bf*beta_hrf(1:size(bf,2),:); %HRF
elseif job.HRFE.hrfm==2 || job.HRFE.hrfm==4
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
    save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'para', 'hrfa', 'event_bold', 'event_number','PARA','-v7.3');
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
    toc
end

if ~flag_ROI % 3D data
    v0 = v(1);   v0.dt = [16,0]; v0.n = [1 1];
    HRF_para_str = {'Height.nii', 'Time2peak.nii','FWHM.nii'};
    dat= zeros(v0.dim);
    if job.savedata.hrfnii_save
        %write nifti maps of the three HRF parameters
        for i=1:3
            v0.fname = fullfile(outdir,[job.prefix,name,'_',HRF_para_str{i}]);
            dat(voxel_ind)=PARA(i,:);
            spm_write_vol(v0,dat);
        end
        % write number of events
        v0.fname = fullfile(outdir,[job.prefix,name,'_event_number.nii']);
        dat(voxel_ind)=event_number;
        spm_write_vol(v0,dat);
    end
    if job.savedata.hrfmat_save
        save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'voxel_ind','v0','-append')
    end
    
    if job.rmoutlier %only for 3D data
        PARA_rm = PARA;
        pvalue_rm=job.para_global.pvalue_rm;
        Inpainted_method =job.para_global.Inpainted ;
        for i=1:3
            v0.fname = fullfile(outdir,[job.prefix,name,'_Olrm_',HRF_para_str{i}]);
            dat(voxel_ind)=PARA(i,:);
            if i==1
                id = (PARA_rm(1,:)<=0);
                PARA_rm(1,id)=nan;
                [PARA_rm(1,:),idx,outliers] = deleteoutliers(PARA_rm(1,:),pvalue_rm,1);
                fprintf('Response height outlier: negative/zeros #%d + extreme #%d, [%3.3f %3.3f];\n',nnz(id),length(idx),min(outliers),max(outliers));
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
            if job.savedata.hrfnii_save
                spm_write_vol(v0,dat);
            end
        end
        if job.savedata.hrfnii_save || job.savedata.hrfmat_save
            % write number of events
            v0.fname = fullfile(outdir,[job.prefix,name,'_Olrm_event_number.nii']);
            dat(voxel_ind)=event_number;
            event_number_rm  = event_number;
            event_number_rm(id_rm)=nan;
            dat(voxel_ind)=event_number_rm;
            dat=inpaint_nans3(dat,Inpainted_method);
            event_number_rm = dat(voxel_ind);
            if job.savedata.hrfnii_save
                spm_write_vol(v0,dat);
            end
            if job.savedata.hrfmat_save
                save(fullfile(outdir,[job.prefix,name, '_hrf.mat']), 'PARA_rm','event_number_rm','id_rm','-append')
            end
        end
    end
    
    if ~isempty(connroinfo) && any(conndata_flag~=2)
        conid = [conndata_flag~=2];
        data_3D = nan(Nscans, prod(v0.dim)); 
        data_3D(:,voxel_ind) = data;
        wgr_conn_run(data_3D, connroinfo(conid,:),v0,name,outdir,flag_pval_pwgc);    
        clear data_3D data
    end


    if job.savedata.deconv_save || any(conndata_flag~=1)
        % writing back deconvolved data into nifti file
        v1 = v; dat3 = zeros(v1(1).dim);
        for i=1:Nscans
            v1(i).fname = fullfile(outdir,[job.prefix,name,'.nii']);
            v1(i).dt = [16,0]; 
            dat3(voxel_ind) = data_deconv(i,:);
            if job.savedata.deconv_save
                spm_write_vol(v1(i),dat3);
            end
        end

        if ~isempty(connroinfo) && any(conndata_flag~=1)
            conid =  [conndata_flag~=1];
            data_deconv3D = nan(Nscans, prod(v0.dim)); 
            data_deconv3D(:,voxel_ind) = data_deconv;
            name2 = [name,'_deconv'];
            wgr_conn_run(data_deconv3D, connroinfo(conid,:),v0,name2,outdir,flag_pval_pwgc);                            
        end

        if job.rmoutlier && job.para_global.rmoutlier_deconv
            v1 = v; dat3 = zeros(v1(1).dim);
            for i=1:Nscans
                v1(i).fname = fullfile(outdir,[job.prefix,name,'_Olrm.nii']);
                v1(i).dt = [16,0]; 
                data_deconv_rm= data_deconv(i,:);
                data_deconv_rm(id_rm)=nan;
                dat3(voxel_ind) = data_deconv_rm;
                dat3=inpaint_nans3(dat3,Inpainted_method);
                data_deconv3D(i,:) = dat3(:);
                if job.savedata.deconv_save
                    spm_write_vol(v1(i),dat3);
                end
            end

            if ~isempty(connroinfo)&& any(conndata_flag~=1)
                conid =  [conndata_flag~=1];
                name2 = [name,'_deconv_Olrm'];
                wgr_conn_run(data_deconv3D, connroinfo(conid,:),v0,name2,outdir,flag_pval_pwgc);     
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
        wgr_deconv_job(job,data,flag_ROI, outdir)
    else
        wgr_conn_job(job,data,flag_ROI, outdir)
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
        [mat,atlas,nii]= wgr_check_ROI(genericROI, job.ref_nii,job.para_global.delete_files);  
        connroinfo{j,1} = [mat;atlas;nii];
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

function wgr_conn_job(job,data,flag_ROI, outdir, v, voxel_ind)
Nscans = size(data,1);
name = job.raw_outname;
flag_pval_pwgc = job.para_global.flag_pval_pwgc;
if ~flag_ROI % 3D data
    v0 = v(1);   v0.dt = [16,0]; v0.n = [1 1];    
    if ~isempty(job.connectivity)
        fprintf('Connectivity analysis...\n');
        connroinfo = wgr_conn_check(job);
    end
    if ~isempty(connroinfo) 
        data_3D = nan(Nscans, prod(v0.dim)); 
        data_3D(:,voxel_ind) = data;
        wgr_conn_run(data_3D, connroinfo,v0,name,outdir,flag_pval_pwgc);    
        clear data_3D data
    else
        warning('Please Configure Parameters for Connectivity Analysis')
    end
else %ROI-wise
    if ~isempty(job.connectivity)
        fprintf('Connectivity analysis...\n');
        connroinfo = wgr_conn_check(job);
        wgr_conn_run(data, connroinfo,[],name,outdir);
        clear data
    else
        warning('Please Configure Parameters for Connectivity Analysis')
    end
end

if job.job_save
    save(fullfile(outdir,[name, '_conn_job.mat']), 'job');  
end

function wgr_conn_run(data, connroinfo,v0,name,outdir,flag_pval_pwgc)
%data: nobs x nvar (3D index)
fprintf('Connectivity analysis...\n ')
meastr = {'pwGC','CGC','Pearson','Spearman','PCGC'};
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
    if nroi==0 && isempty(v0)
        ROI = data;
    end
    if ~flag_seedROI % seed map
        voxel_ind = find(var(data)>0);
        dat.data = data(:,voxel_ind);
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

            [M] = wgr_CGC_OLS(dat,order,ndinfo,conn_type, flag_pval_pwgc);
            ordeinfo = ['_order',num2str(order)];
            for i=1:nroi
                if nroi>1
                    tmp = [num2str(i),'_'];
                else
                    tmp = '';
                end
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_pwGC',ordeinfo,'.nii']);
                gc = M.GC_Matrix_Out(i,:);
                dat3(voxel_ind) = gc;
                spm_write_vol(v0,dat3);
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_N_pwGC',ordeinfo,'.nii']);
                dat3(voxel_ind) = wgr_pwgc_2normal(gc,M.nobs,order);
                spm_write_vol(v0,dat3);
                if flag_pval_pwgc
                    v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_outflow_pval_pwGC',ordeinfo,'.nii']);
                    dat3(voxel_ind) = M.pval_GC_Matrix_Out(i,:);
                    spm_write_vol(v0,dat3);
                end
                
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_pwGC',ordeinfo,'.nii']);
                gc = M.GC_Matrix_In(i,:);
                dat3(voxel_ind) = gc;
                spm_write_vol(v0,dat3);             
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_N_pwGC',ordeinfo,'.nii']);
                dat3(voxel_ind) = wgr_pwgc_2normal(gc,M.nobs,order);
                spm_write_vol(v0,dat3);
                if flag_pval_pwgc
                    v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_inflow_pval_pwGC',ordeinfo,'.nii']);
                    dat3(voxel_ind) = M.pval_GC_Matrix_In(i,:);
                    spm_write_vol(v0,dat3);
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
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_corr_',meastr{connroinfo{j,3}},'.nii']);
                dat3(voxel_ind) = M.Matrix_r(i,:);
                spm_write_vol(v0,dat3);
                v0.fname = fullfile(outdir,[connroinfo{j,6},tmp,name,'_Z_',meastr{connroinfo{j,3}},'.nii']);
                dat3(voxel_ind) = M.Matrix_z(i,:);
                spm_write_vol(v0,dat3);
            end
        end           
    else %ROI to ROI
        dat.data = ROI;
        dat.ROI = [];
        if conn_type==1 || conn_type==2 || conn_type==5 % 1:pairwise  2:conditional 5: partically conditioned
            M = wgr_CGC_OLS(dat,order,ndinfo, conn_type);
            save(fullfile(outdir,[connroinfo{j,6},name,'_',meastr{connroinfo{j,3}},'.mat']),'M');
        else
            M = wgr_FC(dat,conn_type);
            save(fullfile(outdir,[connroinfo{j,6},name,'_Corr_',meastr{connroinfo{j,3}},'.mat']),'M');
        end
        
    end
end

function c = wgr_pwgc_2normal(gc,nobs,order)
c = (nobs-order).*gc - (order-1)/3;
c(c<0)=0;
c = sqrt(c);

function [M] = wgr_CGC_OLS(dat,order,ndinfo,flag_pw_cgc,flag_pval_pwgc, m);
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

if nargin<9
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
                    wgr_seedGC_OLS(ind_X,ind_Y,ROI,data,order,1); 
            else
                [matrix_out{i,j}, matrix_in{i,j}] = wgr_seedGC_OLS(ind_X,ind_Y,ROI,data,order,0); %only pairwise.
            end
        end
    end
else % data to data
    if flag_pw_cgc==5 
        [y_inform_gain, cind] = wgr_init_partial_conditioning(data,[],ndmax,order);
        M.information_gain = y_inform_gain;
        M.condition_id = cind;
    end
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
            if flag_pw_cgc~=5 
                cind=[];
                nd=[];
            end
            [matrix{i,j},p{i,j}] = wgr_CGC_OLS_job(ind_X,ind_Y,data,nvar,order,cind,nd,flag_pw_cgc);
        end
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
    M.GC_Matrix_Out(M.GC_Matrix_Out<0)=0;
    M.GC_Matrix_In(M.GC_Matrix_In<0)=0;
else
    M.GC_Matrix = nan(nvar,nvar);
    for i=1:nbin
        for j=1:nbin
            M.GC_Matrix(indX{i},indY{j}) =  matrix{i,j};
            M.pval_Matrix(indX{i},indY{j}) =  p{i,j};
        end
    end 
    M.GC_Matrix(1:nvar+1:end)=0; %remove diag value.
    M.GC_Matrix(M.GC_Matrix<0)=0;
    if flag_pw_cgc==1
        M.GC_Matrix_N = wgr_pwgc_2normal(M.GC_Matrix,nobs,order);
    end        
end
clear matrix* p_* indX indY

function [gc_out, gc_in, p_out, p_in] = wgr_seedGC_OLS(ind_X,ind_Y,ROI,data1,order,flag_pval);

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
        if flag_pval
            [gc_out(drive,target), p_out(drive,target) ]= wgr_GCA_regress(ROI(:,ind_X(drive)),data1(:,ind_Y(target)), [],order); 
            [gc_in(drive,target), p_in(drive,target)] = wgr_GCA_regress(data1(:,ind_Y(target)),ROI(:,ind_X(drive)), [],order); %Transpose
        else
            gc_out(drive,target) = wgr_GCA_regress(ROI(:,ind_X(drive)),data1(:,ind_Y(target)), [],order); 
            gc_in(drive,target) = wgr_GCA_regress(data1(:,ind_Y(target)),ROI(:,ind_X(drive)), [],order)'; %Transpose
        end
    end
end

function [cgc,pval] = wgr_CGC_OLS_job(ind_X,ind_Y,data1,nvar,order,cind, nd, flag_pw_cgc);
[nvar1]=length(ind_X);
[nvar2]=length(ind_Y);
if flag_pw_cgc==2
    cind = 1:nvar;
end
cgc=zeros(nvar1,nvar2);
pval = nan(nvar1,nvar2);
for drive=1:nvar1
    for target=1:nvar2
        if ind_X(drive) ~= ind_Y(target)
            if flag_pw_cgc==2
                idz = setdiff(cind, [ind_X(drive),ind_Y(target)]);
                conz = data1(:,idz);
            elseif flag_pw_cgc==5
                A = cind(ind_X(drive),:);
                idz = A(~ismembc(A,ind_Y(target)) );
                conz = data1(:,idz(1:nd));
            else
                conz = [];
            end
            [cgc(drive,target),pval(drive,target)] = wgr_GCA_regress(data1(:,ind_X(drive)),data1(:,ind_Y(target)), conz,order); 
        end
    end
end

function [Fy2xIz, pvalue] = wgr_GCA_regress(y,x,z,order);
%% x: N*nx, y: N*ny
%% Fy2x: y->x;
[N, ny] = size(y); 
nz = size(z,2);
%now
X = x(order+1:end,:);

%past
past_ind = repmat([1:order],N-order,1) + repmat([0:N-order-1]',1,order);
X_past = reshape(x(past_ind,:),N-order,order*size(x,2));
Y_past = reshape(y(past_ind,:),N-order,order*size(y,2));
if nz
    Z_past = reshape(z(past_ind,:),N-order,order*size(z,2));
else
    Z_past = [];
end
XZ_past = [X_past  Z_past];
XZY_past = [XZ_past Y_past];

% res_XZ = X - XZ_past*(XZ_past\X);
res_XZ = wgr_regress_res(X,XZ_past);
% res_XZY = X - XZY_past*(XZY_past\X);
res_XZY = wgr_regress_res(X,XZY_past);

%Covariance matrix
cov_XZY= cov(res_XZY);
cov_XZ = cov(res_XZ);

R = det(cov_XZ)/det(cov_XZY);
Fy2xIz = log(R);
if nargout>1
    %nx==1
    m = N-order;      % effective number of observations 
    d1 = order*ny;          % #{full model parameters} - #{reduced model parameters}
    d2 = m-order*(1+ny+nz);       % #{observations} - #{full model parameters}
    mm = d2/d1;
    F = (R-1)*mm;     % F = (exp(gc)-1)*mm;  (RSS_reduced - RSS_full) / RSS_full
    pvalue = 1 - fcdf(F,d1,d2);
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

function resid = wgr_regress_res(y,X)
% copy from function [b,bint,r,rint,stats] = regress(y,X,alpha)
[n,ncolX] = size(X);
% Use the rank-revealing QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);
if isempty(R)
    p = 0;
elseif isvector(R)
    p = double(abs(R(1))>0);
else
    p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
end
if p < ncolX
%     warning(message('stats:regress:RankDefDesignMat'));
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(ncolX,1);
b(perm) = R \ (Q'*y);
yhat = X*b;                     % Predicted responses at each data point.
resid = y-yhat;                     % Residuals.

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