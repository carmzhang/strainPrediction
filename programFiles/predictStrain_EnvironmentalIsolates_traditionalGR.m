%% 
% The script is for the prediction of genetic identity for the library
% of environmental isolates. Here, we do 3 fold cross validation where the 
% full dataset consists of 143 strains has 12 replicates. The validation
% set consists of 3 replicates per strain. 
% The plot for all growth dynamics are also included such that each subplot
% represents all replicates for one strain. Supplementary table 1.11. 
% This script uses growth rate as the features.
% 
%%

clear;
close all;
 
data = importEnvironmentalIsolates_traditionalGR(1);%import data with new labels

% OUTPUT of interest: optimized - optimal hyperparameters (kernel, cost, kkt, sigma)
%                     max_acc - maximum accuracy of cross validation
%                     corresponding to the optimal hyperparameter
%                      validation_acc - accuracy of validation set
[per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = envSVMOpt(data, 12, 3, 3);

disp('Optimal SVM hyperparameters: ')
disp(optimized)
disp('CV accuracy: ')
disp(max_acc)
disp('Validation set accuracy: ')
disp(validation_acc)
