%% 
% The script is for the prediction of genetic identity for the library
% of clinical isolates. Here, we do 4 fold cross validation where the 
% full dataset consists of 203 strains has 4 replicates. Below are examples
% of how to run the predictions. To choose which growth condition to run 
% the prediction for, go to line 50 and pick an option described in the
% comments. This script uses the derivative of growth curves as the features.
% 
%%

close all;
clear;

%import dataset with growth rates as the features
%% CHOOSE: features to do the prediction
% OPTIONS: 'gr' - growth rate curve
%          'gc' - smoothed growth curve
%          'AUC' - area under the curve of the smoothed growth curve
%          'u' - maximum growth rate
%          'u+AUC' - both area under the curve and maximum growth rate
data = importClinicalIsolates('gr');
datasize = size(data);
numFeatures = datasize(2) - 1; %first column is class label
per_condition = numFeatures/4;

%indexes of columns corresponding to growth conditions
one = 2:(2+per_condition-1);%10,000x - 2:99 for gr
two = (2+per_condition):(2+(2*per_condition)-1);%phage - 100:197 for gr
three = (2+(2*per_condition)):(2+(3*per_condition)-1);%100x - 198:295 for gr
four = (2+(3*per_condition)):(2+(4*per_condition)-1);%carbenicillin - 296:393 for gr

%cell: each element in cell is the dataset for one set of growth conditions
data_all = {data(:, one);
            data(:, two);
            data(:, three);
            data(:, four);
            [data(:, one), data(:, two)];
            [data(:, one), data(:, three)];
            [data(:, one), data(:, four)];
            [data(:, two), data(:, three)];
            [data(:, two), data(:, four)];
            [data(:, three), data(:, four)];
            [data(:, one), data(:, two), data(:, three)];
            [data(:, one), data(:, two), data(:, four)];
            [data(:, one), data(:, three), data(:, four)];
            [data(:, two), data(:, three), data(:, four)];
            [data(:, one), data(:, two), data(:, three), data(:, four)];
            };

%% CHOOSE growth condition to do prediction on
% OPTIONS - 1: 10,000x
%           2: phage
%           3: 100x
%           4: carbenicillin
%           5: 10,000x and phage
%           6: 10,000x and 100x
%           7: 10,000x and carbenicillin
%           8: phage and 100x
%           9: phage and carbenicillin
%           10: 100x and carbenicillin
%           11: 10,000x and phage and 100x
%           12: 10,000x and phage and carbenicillin
%           13: 10,000x and 100x and carbenicillin
%           14: phage and 100x and carbenicillin
%           15: 10,000x and phage and 100x and carbenicillin

%ex: growth condition is 10,000x for all predictions below
data = [data(:, 1), data_all{1}];
%EX: to change the growth cindition to phage (2), replace the line above
%with data = [data(:, 1), data_all{2}];

%example 1: 3 fold cross validation with holdout (supplementary table 1.1 
% column 3) using growth rates as the features
% OUTPUT of interest: optimized - optimal hyperparameters (kernel, cost, kkt, sigma)
%                     max_acc - maximum accuracy of cross validation
%                     corresponding to the optimal hyperparameter
% (with holdout ONLY) validation_acc - accuracy of validation set
repl = 4;
fold = 3;
holdout = 1;%with holdout
topk = 1;
kernelOpt = 0;%full hyperparameter search
[per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);

%example 2: 4 fold cross validation (supplementary table 1.1 column 2)
%using growth rates as the features
repl = 4;
fold = 4;
holdout = 0;
topk = 1;
kernelOpt = 0;
[per_kernel_2, avg_per_kernel_2, max_acc_2, optimized_2, validation_acc_2, train_acc_2] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);

%Ex 1 output
disp('3 fold CV with holdout - 10,000x growth condition')
disp('Optimal SVM hyperparameters: ')
disp(optimized)
disp('CV accuracy: ')
disp(max_acc)
disp('Validation set accuracy: ')
disp(validation_acc)

%Ex 2 output
disp('4 fold CV without holdout - 10,000x growth condition')
disp('Optimal SVM hyperparameters: ')
disp(optimized_2)
disp('CV accuracy: ')
disp(max_acc_2)

