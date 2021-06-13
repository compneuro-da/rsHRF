%test in Matlab 2016a, 2018a, 2020a, 2020b.
clc,clear,close all
file = 's1_suj5_525_preproc.mat';
load(file)

flag_plot_signal = 1;
if flag_plot_signal
    xxx = 1/FS2*[1:length(LFP)];
    figure1 = figure('Color','w','units','normalized','position',[0.4 0.3 0.5 0.4]);
    % Create axes
    axes1 = axes('Parent',figure1,...
        'Position',[0.05 0.3 0.9 0.45]);
    hold(axes1,'on');
    plot1 = plot(xxx,[LFP -OPT]);
    set(plot1(1),'DisplayName','LFP envelope');
    set(plot1(2),'DisplayName',' - IOS','LineWidth',1);

    box(axes1,'on');
    set(axes1,'FontName','Cambria','FontSize',16); 
    if 1
        legend(axes1,'show');
    end
    ylim([-0.35 0.5])
    xlim(xxx([1 end]))
%     return
end

%% HRF estimation
TR = 1; 
FS = 1/TR;
fa = FS*100;
%% downsample data
bold = resample(-OPT,fa,FS2*100);
neur = resample(LFP,fa,FS2*100);

%% para
para.TR= TR;
para.T=1;
para.T0=1;
if TR==0.5
    para.len = 24;
elseif TR==1
    para.len = 26;
elseif TR==2
    para.len = 30;
end
dt = para.TR/para.T;
para.dt = dt;
para.lag =1:fix(5/para.dt); %ref. Fig2. in Pan et al. NeuroImage 179,207-214,
para.thr = 1;
para.localK = 2;
para.flag_detrend = 1;
para.Band = [0.01 0.1];
opt  = impulseestOptions('RegulKernel','none');
para.options = opt;
para.temporal_mask = [];

flag_all=1;
if flag_all
    [hrf] = rsHRF_estimation_impulseest_IO(neur, bold, para);
    rsHRF_plot_hrf_mean_std_4_LFP_IOS(hrf.H_LFP, hrf.H_LFP_bin, hrf.H_Blind,TR);
    rsHRF_plot_hrf_mean_std_4_LFP_IOS(hrf.H_LFP, hrf.H_LFP_bin, hrf.H_Blind,TR,1);
%     return
end

flag_sliding_windows=0;
if flag_sliding_windows
    HRFall = {};
    h=1;
    sw = 150:10:240;
    for i= sw
        nobs = length(neur);
        k = 1;
        block = i*FS;
        index={};
        for mv = 30*FS:fix(4*FS):nobs-block
            index{k}= [mv: mv+block];
            if 0
                disp([mv mv+block])
            end
            k=k+1;
        end
        fprintf('# sliding windows: %d\n',k-1)
        HRF={};
        for ii=1:length(index)
            HRF{ii,1}  = rsHRF_estimation_impulseest_IO(neur(index{ii}), bold(index{ii}), para);
        end
        HRFall{h,1} = HRF;
        
        for x=1:length(HRF)
            hrf_LFP_bin(:,x) = HRF{x}.H_LFP_bin;
            hrf_LFP(:,x) = HRF{x}.H_LFP;
            hrf_Blind(:,x) =HRF{x}.H_Blind; 
        end

        x=1:length(HRF);
        PARA = arrayfun(@(x) rsHRF_get_HRF_parameters(hrf_Blind(:,x),dt),1:size(hrf_Blind,2), 'UniformOutput', false);
        PARA_Blind = cell2mat(PARA);
        PARA = arrayfun(@(x) rsHRF_get_HRF_parameters(hrf_LFP_bin(:,x),dt),1:size(hrf_LFP_bin,2), 'UniformOutput', false);
        PARA_LFP_bin = cell2mat(PARA);
        PARA = arrayfun(@(x) rsHRF_get_HRF_parameters(hrf_LFP(:,x),dt),1:size(hrf_LFP,2), 'UniformOutput', false);
        PARA_LFP = cell2mat(PARA);
        h1 = [PARA_Blind(1,:)' PARA_LFP_bin(1,:)' PARA_LFP(1,:)'];
        fwhm = [PARA_Blind(2,:)' PARA_LFP_bin(2,:)' PARA_LFP(2,:)'];
        t2p = [PARA_Blind(3,:)' PARA_LFP_bin(3,:)' PARA_LFP(3,:)'];

        [rh,ph]=corr(h1);
        [rw,pw]=corr(fwhm);
        [rt,pt]=corr(t2p);
        HRFall{h,2} = {rh,ph; rw,pw; rt,pt};
        HRFcorr(h,:) = [rh(1,2:3) rh(2,3)  ph(1,2:3) ph(2,3)];
        h=h+1;
    end
    if ~isempty(HRFall)
        save(['fa=',num2str(fa),'_len=',num2str(para.len),'_HRF.mat'], 'HRFall', 'HRFcorr', 'sw', 'FS');
    end
end

flag_plot_slidingwindows = 0;
if flag_plot_slidingwindows
    load(['fa=',num2str(fa),'_len=',num2str(para.len),'_HRF.mat'])
    %% Table 1
    disp('Mean & STD')
    mean(HRFcorr)
    std(HRFcorr)
    % return 
    id = 6;
    HRF = HRFall{id};
    fprintf('windows %d r=%1.4f (LFP-pp), %1.4f (LFP)\n',sw(id),HRFcorr(id,1:2))
    hrf_LFP_bin=[]; hrf_LFP=[]; hrf_Blind=[];
    for x=1:length(HRF)
        hrf_LFP_bin(:,x) = HRF{x}.H_LFP_bin;
        hrf_LFP(:,x) = HRF{x}.H_LFP;
        hrf_Blind(:,x) =HRF{x}.H_Blind; 
    end

    rsHRF_plot_hrf_mean_std_4_LFP_IOS(hrf_LFP,hrf_LFP_bin,hrf_Blind,dt)
    rsHRF_plot_hrf_mean_std_4_LFP_IOS(hrf_LFP,hrf_LFP_bin,hrf_Blind,dt,1)
end
