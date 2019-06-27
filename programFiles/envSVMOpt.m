function [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = envSVMOpt(data, repl, fold, holdoutSize)
%FUNCTION: Optimizes hyperparameters of svm with holdout and returns one  
%          of the optimized parameter sets along with the accuracy using
%          this parameter set (training and validation)
% INPUT:
%   data - first column is label, other columns are features
%   repl - number of replicates per strain
%   fold - n corresponding to n-fold cross validation
%   holdoutSize - size of hold out for validation set in
%   terms of number of replicates per strain
% OUTPUT: 
%   per_kernel - cell with the accuracy and predicted labels of the test
%                set for all hyperparameter sets
%   avg_per_kernel - average accuracy  across the n folds for  each 
%                    parameter set, grouped  by kernel with the same order
%                    as per_kernel
%   max_acc - maximum accuracy out of all hyperparameter sets
%   optimized - the optimal hyperparameters
%   validation_acc - the accuracy for the holdout dataset using one of the 
%                    optimized hyperparameter set
%   train_acc - the accuracy for the training dataset using one of the 
%                    optimized hyperparameter set


s = size(data);
idx_all = 1:s(1);
data_val_idx = [];
for i = 1:holdoutSize
    data_val_idx = [data_val_idx, i:repl:s(1)];
end
data_val_idx = sort(data_val_idx);%indexes of holdout
idx_all(data_val_idx) = [];
s_folds = size(idx_all);

idx_folds = cell(fold, 2);%idx for test (1) and train (2) for each fold
for i = 1:1:fold
    start = (i*3) + 1;
    temp = [start:repl:s(1), ...
            (start+1):repl:s(1), ...
            (start+2):repl:s(1)];
    temp = sort(temp);
    temp_2 = setdiff(idx_all, temp);
    idx_folds{i, 1} = temp;
    idx_folds{i, 2} = temp_2;
end

%% optimize parameter selection for svm
kernelOpt = [string('quadratic'), string('rbf')];
num_kernel = 2;
sigmaGradient = [1, 5, 10, 20];%rbf kernel only
costGradient = [10, 100, 1000];
kktviolationlevelGradient = [0, 0.05, 0.1];

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
    per_kernel{k, 2} = cell(p_size, 3);
    
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
            K = s(1);%number of classes
            M = 98;%dimension of features
            
            train_label = data(idx_folds{i, 2}, 1);
            test_label = data(idx_folds{i, 1}, 1);
            test_feature = testX;
            train_feature = trainX;
            
            %%
            [predict_label,result]=multi_class_svm_cv(test_label,train_label,test_feature,train_feature,N_test,N_train, kernel, cost, kkt, sigma);
            % predict_label.svm is the predicted labels based on Matlab default SVM
            % result.svm is the accuracy of svm
            
            m = [m, result.svm];%store accuracy
            saveScores = [saveScores, predict_label.svm'];
            disp(i)
            
        end
        
        per_kernel{k, 2}{j, 2} = m;
        per_kernel{k, 2}{j, 3} = saveScores;
        
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

