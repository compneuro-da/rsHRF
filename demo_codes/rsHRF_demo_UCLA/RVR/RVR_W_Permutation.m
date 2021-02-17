
function RVR_W_Permutation(Data_Path, Scores, Perm_times_Range, RandIndex_Folder, Covariates, Pre_Method, ResultantFolder, Queue)

TaskQuantity = 100;
JobsPerTask = fix(length(Perm_times_Range) / TaskQuantity);
JobsRemain = mod(length(Perm_times_Range), TaskQuantity * JobsPerTask);

mkdir([ResultantFolder filesep 'rand_perm']);

for i = 1:100
    i
    
    ID = Perm_times_Range([(i - 1) * JobsPerTask + 1 : i * JobsPerTask]);
    if i == 100 & JobsRemain
        ID = [ID Perm_times_Range(length(Perm_times_Range) - JobsRemain + 1 : end)];
    end
    for j = 1:length(ID)
%         Rand_Index = randperm(length(Scores));
        load([RandIndex_Folder filesep 'RandID_' num2str(ID(j)) '.mat']);
        Rand_Index = RandID;
        Rand_Score{j} = Scores(Rand_Index);
    end
    
    save([ResultantFolder filesep 'rand_perm' filesep 'rand_perm_' num2str(i) '.mat'], 'Rand_Score');
    Job_Name = ['perm_W_' num2str(i)];
    pipeline.(Job_Name).command = 'W_Calculate_RVR_SGE(opt.para1, opt.para2, opt.para3, opt.para4, opt.para5, opt.para6)';
    pipeline.(Job_Name).opt.para1 = Data_Path;
    pipeline.(Job_Name).opt.para2 = Rand_Score;
    pipeline.(Job_Name).opt.para3 = ID;
    pipeline.(Job_Name).opt.para4 = Covariates;
    pipeline.(Job_Name).opt.para5 = Pre_Method;
    pipeline.(Job_Name).opt.para6 = ResultantFolder;
    clear Rand_Score;
end

Pipeline_opt.mode = 'qsub';
Pipeline_opt.qsub_options = Queue;
Pipeline_opt.mode_pipeline_manager = 'batch';
Pipeline_opt.max_queued = 200;
Pipeline_opt.flag_verbose = 1;
Pipeline_opt.flag_pause = 0;
Pipeline_opt.flag_update = 1;
Pipeline_opt.path_logs = [ResultantFolder filesep 'logs'];

psom_run_pipeline(pipeline,Pipeline_opt);


