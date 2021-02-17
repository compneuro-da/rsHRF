clc,clear
addpath('E:\rsHRF_demo_UCLA\RVR')
matfile = {
'FC20.mat'    
'FC79.mat'    
};
flag_model_generalization = 1;
if flag_model_generalization
    
    load(matfile{1})
    age20=age;
    fcd20=all_matd; 
    fcr20=all_matr; 
    cov20=cova;
    
    load(matfile{2})
    age79=age;
    fcd79=all_matd; 
    fcr79=all_matr; 
    cov79=cova;


    Pre_Method = 'Scale';
    np=1000;
    r = zeros(np+1,2);
    rd = r; MAE = r; MAEd=r;

    seed = 10;
    rng(seed)
    for i=np+1:-1:1
        RandIndex20(i,:) = randperm(20);
        RandIndex79(i,:) = randperm(79);
    end

    for i=1:np+1
        fprintf('%d',i)
        if i==1
            Permutation_Flag=0;
        else
            Permutation_Flag=1;
        end

        disp('BOLD')
        Prediction79to20{i} = RVR_APredictB_modified(fcr79, age79, cov79, fcr20, age20, cov20, Pre_Method, Permutation_Flag,RandIndex79(i,:));
        Prediction20to79{i} = RVR_APredictB_modified(fcr20, age20, cov20, fcr79, age79, cov79, Pre_Method, Permutation_Flag,RandIndex20(i,:));
        r(i,1) = Prediction79to20{i}.Corr; 
        r(i,2) = Prediction20to79{i}.Corr; 
        MAE(i,1) = Prediction79to20{i}.MAE; 
        MAE(i,2) = Prediction20to79{i}.MAE;  

        disp('Deconvolved BOLD')
        dPrediction79to20{i} = RVR_APredictB_modified(fcd79, age79, cov79, fcd20, age20, cov20, Pre_Method, Permutation_Flag,RandIndex79(i,:));
        dPrediction20to79{i} = RVR_APredictB_modified(fcd20, age20, cov20, fcd79, age79, cov79, Pre_Method, Permutation_Flag,RandIndex20(i,:));
        rd(i,1) = dPrediction79to20{i}.Corr; 
        rd(i,2) = dPrediction20to79{i}.Corr; 
        MAEd(i,1) = dPrediction79to20{i}.MAE; 
        MAEd(i,2) = dPrediction20to79{i}.MAE;  
    end

    p_r = sum( r(2:end,:)>=repmat(r(1,:),np,1))./np;
    p_MAE = sum( MAE(2:end,:)<=repmat(MAE(1,:),np,1))./np;

    p_rd = sum( rd(2:end,:)>=repmat(rd(1,:),np,1))./np;
    p_MAEd = sum( MAEd(2:end,:)<=repmat(MAEd(1,:),np,1))./np;

    save(['model_generalization_',num2str(np),'_rng',num2str(seed),'.mat'],'*Predic*','p_*','MAE*','r','rd');
    
    
else % prediction in each dataset.
    
    
    for j=1:2
        load(matfile{j})
        Pre_Method = 'Normalize';
        np=1000;
        seed = 10;
        rng(seed)
        for kkk=1:2
            if kkk==1
                disp('raw FC')
                FCmat = all_matr;
            else
                disp('Deconv FC')
                FCmat = all_matd;
            end
            all_behav = age;
            kf = length(all_behav);
            for i=1:np+1
                fprintf('%d ',i)
                if i==1
                    Prediction = RVR_LOOCV_modified(FCmat,all_behav', cova, Pre_Method, 0, 0);
                else
                    Prediction = RVR_LOOCV_modified(FCmat,all_behav', cova, Pre_Method, 0, 1);
                end
                rall(i,kkk) = Prediction.Corr; 
                MAE(i,kkk) = Prediction.MAE;
                y_predict_all{i,kkk} =Prediction;
            end
        end
        p_r = sum( rall(2:end,:)>=repmat(rall(1,:),np,1))./np;
        p_MAE = sum( MAE(2:end,:)<=repmat(MAE(1,:),np,1))./np;
        save(['FC',num2str(length(age)),'_permutation_',num2str(np),'_rng',num2str(seed),'.mat'],'rall','y_predict_all','p_*','MAE*')
    end
end