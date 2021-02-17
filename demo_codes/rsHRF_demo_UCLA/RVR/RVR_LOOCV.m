function Prediction = RVR_LOOCV(Subjects_Data, Subjects_Scores, Covariates, Pre_Method, Weight_Flag, ResultantFolder)
%
% Subject_Data:
%           m*n matrix
%           m is the number of subjects
%           n is the number of features
%
% Subject_Scores:
%           the continuous variable to be predicted, [1*m]
%
% Covariates:
%           m*n matrix
%           m is the number of subjects
%           n is the number of covariates
%
% Pre_Method:
%          'Normalize', 'Scale', 'None'
%
% Weight_Flag:
%           whether to compute the weight, 1 or 0
%
% ResultantFolder:
%           the path of folder storing resultant files
%

if nargin >= 5
    if ~exist(ResultantFolder, 'dir')
        mkdir(ResultantFolder);
    end
end

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);

for i = 1:Subjects_Quantity
    
    disp(['The ' num2str(i) ' subject!']);
    
    Training_data = Subjects_Data;
    Training_scores = Subjects_Scores;
    
    % Select training data and testing data
    test_data = Training_data(i, :);
    test_score = Training_scores(i);
    Training_data(i, :) = [];
    Training_scores(i) = [];
    
    if ~isempty(Covariates)
        Covariates_test = Covariates(i, :);
        Covariates_training = Covariates;
        Covariates_training(i, :) = [];
        [Training_quantity, Covariates_quantity] = size(Covariates_training);
        M = 1;
        for j = 1:Covariates_quantity
            M = M + term(Covariates_training(:, j));
        end
        slm = SurfStatLinMod(Training_data, M);
        
        Training_data = Training_data - repmat(slm.coef(1, :), Training_quantity, 1);
        for j = 1:Covariates_quantity
            Training_data = Training_data - ...
                repmat(Covariates_training(:, j), 1, Feature_Quantity) .* repmat(slm.coef(j + 1, :), Training_quantity, 1);
        end
    end
    
    if strcmp(Pre_Method, 'Normalize')
        %Normalizing
        MeanValue = mean(Training_data);
        StandardDeviation = std(Training_data);
        [~, columns_quantity] = size(Training_data);
        for j = 1:columns_quantity
            Training_data(:, j) = (Training_data(:, j) - MeanValue(j)) / StandardDeviation(j);
        end
    elseif strcmp(Pre_Method, 'Scale')
        % Scaling to [0 1]
        MinValue = min(Training_data);
        MaxValue = max(Training_data);
        [~, columns_quantity] = size(Training_data);
        for j = 1:columns_quantity
            Training_data(:, j) = (Training_data(:, j) - MinValue(j)) / (MaxValue(j) - MinValue(j));
        end
    end
    Training_data_final = double(Training_data);

    % Covariate test data
    if ~isempty(Covariates)
        test_data = test_data - slm.coef(1, :);
        for j = 1:Covariates_quantity
            test_data = test_data - repmat(Covariates_test(j), 1, Feature_Quantity) .* slm.coef(j + 1, :);
        end
    end
    % Normalize test data
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        test_data = (test_data - MeanValue) ./ StandardDeviation;
    elseif strcmp(Pre_Method, 'Scale')
        % Scale
        test_data = (test_data - MinValue) ./ (MaxValue - MinValue);
    end
    test_data_final = double(test_data);

    % RVR training & predicting
    d.train{1} = Training_data_final * Training_data_final';
    d.test{1} = test_data_final * Training_data_final';
    d.tr_targets = Training_scores';
    d.use_kernel = 1;
    d.pred_type = 'regression';
    output = prt_machine_rvr(d, []);
    
    Predicted_Scores(i) = output.predictions; 
    
end

Prediction.Score = Predicted_Scores;
[Prediction.Corr, ~] = corr(Predicted_Scores', Subjects_Scores');
Prediction.MAE = mean(abs((Predicted_Scores - Subjects_Scores)));

if nargin >= 5
    save([ResultantFolder filesep 'Prediction_res.mat'], 'Prediction');
    disp(['The correlation is ' num2str(Prediction.Corr)]);
    disp(['The MAE is ' num2str(Prediction.MAE)]);
    % Calculating w
    if Weight_Flag
        W_Calculate_RVR(Subjects_Data, Subjects_Scores, Covariates, Pre_Method, ResultantFolder); 
    end
end