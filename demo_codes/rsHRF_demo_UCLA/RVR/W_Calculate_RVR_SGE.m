
function W_Calculate_RVR_SGE(Subjects_Data_Path, Rand_Scores, ID, Covariates, Pre_Method, ResultantFolder)

tmp = load(Subjects_Data_Path);
FieldName = fieldnames(tmp);

for i = 1:length(ID)
    disp(ID(i));
    [w_Brain, ~] = W_Calculate_RVR(tmp.(FieldName{1}), Rand_Scores{i}, Covariates, Pre_Method);
    save([ResultantFolder filesep 'w_Brain_' num2str(ID(i)) '.mat'], 'w_Brain');
end