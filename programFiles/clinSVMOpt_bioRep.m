function [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt_bioRep(data, repl, fold, holdout, topk, kOpt)
%FUNCTION: Optimizes hyperparameters of svm with holdout and returns one
%          of the optimimal parameter sets. If holdout is done then the
%          accuracy of the validation and training set are returned using
%          one of the optimal parameter sets.
% INPUT:
%   data - first column is label, other columns are features
%   repl - number of replicates per strain
%   fold - n corresponding to n-fold cross validation (ex: n = 3 for 3
%   fold cross validation)
%       OPTIONS - 2 or 3
%   holdout - hold out for validation set (0 for no holdout) 1 for holdout
%       OPTIONS - 0 or 1
%   topk - top k predictions are looked at for correct class
%   kOpt - kernel options for hyperparameter search
%       OPTIONS - 0 for linear, quadratic, and rbf
%               - 1 for rbf only
% OUTPUT:
%   per_kernel - cell with the accuracy and predicted labels of the test
%                set for all hyperparameter sets
%   avg_per_kernel - average accuracy  across the n folds for  each
%                    parameter set, grouped  by kernel with the same order
%                    as per_kernel
%   max_acc - maximum accuracy out of all hyperparameter sets
%   optimized - the optimal hyperparameters
%   validation_acc - the accuracy for the holdout dataset using one of the
%                    optimized hyperparameter set (only for holdout)
%   train_acc - the accuracy for the training dataset using one of the
%                    optimized hyperparameter set (only for holdout)

s = size(data);
idx_all = 1:s(1);

if(holdout == 1)
    data_val_idx = [];
    data_val_idx = [[9:repl:s(1)], [10:repl:s(1)], [11:repl:s(1)], [12:repl:s(1)]];
    data_val_idx = sort(data_val_idx);%indexes of holdout
    idx_all(data_val_idx) = [];
    
    idx_folds = cell(fold, 2);%idx for test (1) and train (2) for each fold
    for i = 1:1:fold
        start = (i-1)*4 + 1;        
        temp = [[start:repl:s(1)], [start+1:repl:s(1)], [start+2:repl:s(1)], [start+3:repl:s(1)]];
        temp = sort(temp);
        temp_2 = setdiff(idx_all, temp);
        idx_folds{i, 1} = temp;
        idx_folds{i, 2} = temp_2;
    end
elseif(holdout == 0)
    
    idx_folds = cell(fold, 2);%idx for test (1) and train (2) for each fold
    for i = 1:1:fold
        start = (i-1)*4 + 1;        
        temp = [[start:repl:s(1)], [start+1:repl:s(1)], [start+2:repl:s(1)], [start+3:repl:s(1)]];
        temp = sort(temp);
        temp_2 = setdiff(idx_all, temp);
        idx_folds{i, 1} = temp;
        idx_folds{i, 2} = temp_2;
    end
end

%% optimize parameter selection for svm
if(kOpt == 0)
    kernelOpt = [string('linear'), string('quadratic'), string('rbf')];
    num_kernel = 3;
    sigmaGradient = [1, 10, 20];%rbf kernel only
    costGradient = [10, 100, 1000];
    kktviolationlevelGradient = [0, 0.05, 0.1];
elseif(kOpt == 1)
    kernelOpt = [string('rbf')];
    num_kernel = 1;
    sigmaGradient = [1, 10, 20];%rbf kernel only
    costGradient = [10, 100];
    kktviolationlevelGradient = [0, 0.05, 0.1];
elseif(kOpt == 2)
    kernelOpt = [string('rbf')];
    num_kernel = 1;
    sigmaGradient = [10, 20];%rbf kernel only
    costGradient = [10, 100];
    kktviolationlevelGradient = [0.1];
end

%% prediction
per_kernel = cell(num_kernel, 2);%store accuracy and predictions
for k = 1:1:num_kernel
    kernel = kernelOpt(k);
    %combination of parameters
    if(kernel ~= 'rbf') %1. for linear and quadratic kernels
        paramAll = combvec(costGradient, kktviolationlevelGradient);
    else %2. for rbf kernel
        paramAll = combvec(costGradient, kktviolationlevelGradient, sigmaGradient);
    end
    
    p_size = length(paramAll);
    
    per_kernel{k, 1} = kernel;
    per_kernel{k, 2} = cell(p_size, 2);
    
    for j = 1:1:p_size
        
        m = [];
        saveScores = [];
        
        curr_param = paramAll(:,j);
        cost = curr_param(1);
        kkt = curr_param(2);
        
        %store parameter combination
        if(kernel ~= 'rbf')
            sigma = 0;%no sigma option
            param = strcat('kernel: ', kernel, ', cost: ', num2str(cost), ', kkt: ', num2str(kkt));
        else
            sigma = curr_param(3);
            param = strcat('kernel: ', kernel, ', cost: ', num2str(cost), ', kkt: ', num2str(kkt), ', sigma: ', num2str(sigma));
        end
        
        per_kernel{k, 2}{j, 1} = param;
        
        for i = 1:1:fold
            
            trainX = data(idx_folds{i, 2}, 2:end);
            testX = data(idx_folds{i, 1}, 2:end);
            
            N_train = length(idx_folds{1, 2});% number of training samples
            N_test = length(idx_folds{1, 1});% number of testing samples
            
            train_label = data(idx_folds{i, 2}, 1);
            test_label = data(idx_folds{i, 1}, 1);            
            test_feature = testX;
            train_feature = trainX;
            
            %%
            [result]=multi_class_svm_cv_topk(test_label,train_label,test_feature,train_feature,N_test,N_train, kernel, cost, kkt, sigma, topk);
            
            m = [m, result.svm];%store accuracy
            disp(i)
            
        end
        
        per_kernel{k, 2}{j, 2} = m;
        
        disp('...')
        disp(j)
        
    end %end params
    
    disp('... ...')
    disp(k)
    
end%end kernel options

%take average across folds and find optimal parameter set
avg_per_kernel = cell(num_kernel, 2);
for i = 1:1:num_kernel
    avg_per_kernel{i, 1} = per_kernel{i, 1};
    s = size(per_kernel{i, 2});
    all = [];
    for j = 1:1:s(1)
        all = [all, mean(per_kernel{i, 2}{j, 2})];
    end
    avg_per_kernel{i, 2} = all;
end

max_per_kernel = [];
for i = 1:num_kernel
    max_per_kernel = [max_per_kernel; max(avg_per_kernel{i, 2})];
end
max_acc = max(max_per_kernel);

%find best hyperparameters
optimized = [];
for i = 1:1:num_kernel
    idx_max = find(avg_per_kernel{i, 2} == max_acc);
    if(~isempty(idx_max))
        optimized = [optimized; per_kernel{i, 2}(idx_max, 1)];
    end
end

validation_acc = -1;
train_acc = -1;

if(holdout == 1)
    %test optimized hyperparameters on validation set
    validation = data(data_val_idx, 2:end);
    validation_label = data(data_val_idx, 1);
    
    train = data(idx_all, 2:end);
    train_label = data(idx_all, 1);
    
    N_test = length(data_val_idx);
    N_train = length(idx_all);
    
    kernel = extractBetween(optimized{1},':', ',');
    kernel = kernel(1);
    param_num = regexp(optimized{1},'\d*','Match');
    if(strcmp('rbf', kernel))
        cost = str2double(param_num(1));
        kkt = str2double(param_num(2));
        sigma = str2double(param_num(3));
    else
        cost = str2double(param_num(1));
        kkt = str2double(param_num(2));
        sigma = 0;%no parameter for non-rbf kernels
    end
    [out,predict_label,result, F, out2, predict_label2, result2] = multi_class_svm_ci(validation_label, train_label, validation, train, N_test, N_train, kernel, cost, kkt, sigma);
    validation_acc = result.svm;
    train_acc = result2.svm;
end

end

