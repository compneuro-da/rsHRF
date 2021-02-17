clc,clear,close all
addpath E:\data
matdir= 'E:\data\UCLA_ROI9_FC_out';
mo={'raw FC','Deconv FC'};
pvalue={};
flag_run_pvalue = 0; % it may take several hours.
if 1 %scatter plot & p-value
    for scanner = 1:2
        disp(['scanner ',num2str(scanner)])
        if scanner==1
            load(fullfile(matdir,'scanner1_pmask.mat'))
        elseif scanner==2
            load(fullfile(matdir,'scanner2_pmask.mat'))
        end
        tit={'Positive','Negative'};%'All',
        for kkk=1:2
            disp(mo{kkk})
            r=rall{kkk};
            y_predict0 = y_predict_all{kkk};
            pp = [0.001:0.001:0.5];
            for i=1:2
                disp(tit{i})
                rr = r(:,i); 
                [sr,order] = sort(rr,'descend');
                pthresh = pp(order(1));
                flag_pvalue=1;
                if flag_pvalue
                    all_matr = datz;
                    all_matd = datzd;
                    all_behav = age;
                    nsub = length(all_behav);
                    hhh=1;performance_rand=[];ord=[];
                    kf = nsub;
                    if flag_run_pvalue
                        while hhh<5001
                            ord(hhh,:) = randperm(nsub);
                            all_behav0 = all_behav(ord(hhh,:));
                            try
                                [y_predictr, performance0, pmask_kfoldr] = rsHRF_CPM(all_matr,all_behav0,kf,pthresh,cova);
                                performance_rand(:,hhh) = performance0(i,1);
                                hhh=hhh+1;
                            end
                        end
                        randorder{scanner,kkk,i} = ord;
                        performance_randraw{scanner,kkk,i} = performance_rand;
                        pvalue{scanner,kkk,i} = sum(performance_rand>sr(1), 2)./size(performance_rand,2);
                    end
                end
            end
        end   
    end
    if flag_run_pvalue
        save(fullfile(matdir,'pvalue_CPM.mat'),'pvalue','randorder','performance_randraw')
    end
end

if 1 % model generalization
    for scanner=1:2
        if scanner==1
            disp('predict scannner2 from scannner1 FC features')
            load(fullfile(matdir,'scanner1_pmask.mat'))
            pm_raw0 = pm_raw;
            pm_deconv0 = pm_deconv;
            load(fullfile(matdir,'scanner2_pmask.mat'))
            all_matr = datz;
            all_matd = datzd;
            all_behav = age;
        else
            disp('predict scannner1 from scannner2 FC features')
            load(fullfile(matdir,'scanner2_pmask.mat'))
            pm_raw0 = pm_raw;
            pm_deconv0 = pm_deconv;
            load(fullfile(matdir,'scanner1_pmask.mat'))
            all_matr = datz;
            all_matd = datzd;
            all_behav = age;
        end
        nsub = length(all_behav);
        kf = nsub;
        pthresh=2;
        tit={'Positive','Negative'};
        for i=1:2
            disp(tit{i})
            [y_predictr, performancer, pmask_kfoldr] = rsHRF_CPM(all_matr,all_behav,kf,pthresh,cova,pm_raw0{i});
            [y_predictd, performanced, pmask_kfoldd] = rsHRF_CPM(all_matd,all_behav,kf,pthresh,cova,pm_deconv0{i});
            disp(mo{1})
            performance_raw(scanner,i)=performancer(i,1);
            disp(mo{2})
            performance_deconv(scanner,i)=performanced(i,1);
            
            if 1
                figure;
                subplot(1,2,1)
                scatter(age,y_predictr(:,i))
                xlabel('Actual Age'); ylabel('Predicted Age')
                rr(i)=corr(age,y_predictr(:,i));
                title(sprintf('BOLD (r=%1.2f)',rr(i)))
                subplot(1,2,2)
                scatter(age,y_predictd(:,i));
                xlabel('Actual Age'); ylabel('Predicted Age')
                rd(i) = corr(age,y_predictd(:,i));
                title(sprintf('Deconvolved BOLD (r=%1.2f)',rd(i)))
               suptitle(tit{i})
            end
            if flag_run_pvalue
                hhh=1;performance_randr=[];ord=[];
                while hhh<5001
                    ord(hhh,:) = randperm(nsub);
                    all_behav0 = all_behav(ord(hhh,:));
                    try
                        [y_predictr, performance0, pmask_kfoldr] = rsHRF_CPM(all_matr,all_behav0,kf,pthresh,cova,pm_raw0{i});
                        performance_randr(:,hhh) = performance0(i,1);
                        hhh=hhh+1;
                    end
                end
                randorderraw{scanner,i} = ord;
                performance_randraw{scanner,i} = performance_randr;
                pvaluer{scanner,i} = sum(performance_randr>performancer(i,1), 2)./size(performance_randr,2);

                hhh=1;performance_randd=[];ord=[];
                while hhh<5001
                    ord(hhh,:) = randperm(nsub);
                    all_behav0 = all_behav(ord(hhh,:));
                    try
                        [y_predictr, performance0, pmask_kfoldr] = rsHRF_CPM(all_matd,all_behav0,kf,pthresh,cova,pm_deconv0{i});
                        performance_randd(:,hhh) = performance0(i,1);
                        hhh=hhh+1;
                    end
                end
                randorderdeconv{scanner,i} = ord;
                performance_randdeconv{scanner,i} = performance_randd;
                pvalued{scanner,i} = sum(performance_randd>performanced(i,1), 2)./size(performance_randd,2);
            end
        end
    end
    save(fullfile(matdir,'model_gene_CPM.mat'),'performance_raw','performance_deconv');
    if flag_run_pvalue
        save(fullfile(matdir,'pvalue_gene_CPM.mat'),'pvaluer','pvalued','performance_randdeconv',...
            'performance_randraw','randorderdeconv','randorderraw');
    end
end