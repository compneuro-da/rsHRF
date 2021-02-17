function rsHRF_plot_hrf_mean_std_4_LFP_IOS(hrfraw,hrfpp,hrfe,dt,flag_norm)
if nargin<5
    flag_norm=0;
end
% fig1 = figure('color','w');
nobs = size(hrfraw,1);
xx = (1:nobs)*dt;
figure1 = figure('color','w');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% LFP
mf0 = mean(hrfraw,2);
if flag_norm
    mf0 = mf0./max(mf0);
end
sf0 = std(hrfraw,0,2);
plot(xx,mf0,'LineWidth',2,'color','g');%, 'LineSmoothing','on'
%LFP-pp
mf = mean(hrfpp,2);
sf = std(hrfpp,0,2); 
if flag_norm
    mf = mf./max(mf);
end
plot(xx,mf,'LineWidth',2,'color','b');%, 'LineSmoothing','on'

%blind HRF
mf2 = mean(hrfe,2);
sf2 = std(hrfe,0,2); 
if flag_norm
    mf2 = mf2./max(mf2);
end
plot(xx,mf2,'LineWidth',1,'color','r');%, 'LineSmoothing','on'
if any(sf0)
    %LFP
    fill( [xx fliplr(xx)],  [mf0(:,1)'+sf0(:,1)' fliplr(mf0(:,1)'-sf0(:,1)')],'g', 'LineStyle','none','LineWidth',1,'EdgeColor','g'); %, 'LineStyle','-.'
    %LFP-pp
    fill( [xx fliplr(xx)],  [mf(:,1)'+sf(:,1)' fliplr(mf(:,1)'-sf(:,1)')],'b', 'LineStyle','none','LineWidth',1,'EdgeColor','b');
    %blind HRF
    fill( [xx fliplr(xx)],  [mf2(:,1)'+sf2(:,1)' fliplr(mf2(:,1)'-sf2(:,1)')],'r', 'LineStyle','none','LineWidth',1,'EdgeColor','r');
    alpha(0.3)
end
xlim([0 max(xx)])
axis square
% ylim([min([min(mf0-sf0) min(mf-sf) min(mf2-sf2)])*0.95,max([max(mf0+sf0) max(mf+sf) max(mf2+sf2)])*1.1])
ylim([min([min(mf0-sf0) min(mf-sf) min(mf2-sf2)])*1.1,max([max(mf0+sf0) max(mf+sf) max(mf2+sf2)])*1.1])

legend({'HRF_{LFP}','HRF_{LFP-pp}','HRF_{Blind}'},'box','off')
set(axes1,'FontName','Cambria','FontSize',20); %'Calibri','Times New Roman'
box on