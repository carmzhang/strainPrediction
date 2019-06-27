%%
% The script is for the prediction of genetic identity for the library
% of clinical isolates. Here, we do 4 fold cross validation where the
% full dataset consists of 203 strains has 4 replicates. We change the time
% span and plot how this affects the predictive accuracy. Supplementary
% Figure 1.1. This script uses the derivative of growth curves as the features.
%
%%

close all;
clear;

%import dataset with growth rates as the features
data = importClinicalIsolates('gr');
datasize = size(data);
numFeatures = datasize(2) - 1; %first column is class label
per_condition = numFeatures/4;

%time range to compare
time_range = [98, 88, 78, 68, 58, 48];

data_all = cell(15, length(time_range));
for i = 1:length(time_range)
    
    per_condition = time_range(i);
    %indexes of columns corresponding to growth conditions
    one = 2:(2+per_condition-1);%10,000x - 2:99 for gr
    two = 100:(100+per_condition-1);%phage - 100:197 for gr
    three = 198:(198+per_condition-1);%100x - 198:295 for gr
    four = 296:(296+per_condition-1);%carbenicillin - 296:393 for gr
    
    %cell: each element in cell is the dataset for one set of growth conditions
    data_all(:, i) = {data(:, one);
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
end

%% growth conditions to do prediction on
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

%4 fold cross validation using growth rates as the features
repl = 4;
fold = 4;
holdout = 0;
topk = 1;
kernelOpt = 1;

%% Single growth conditions
all_acc = zeros(4, length(time_range));
for i = 1:4
    for j = 1:length(time_range)
        data = [data(:, 1), data_all{i, j}];
        [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);
        all_acc(i, j) = max_acc;
    end
end

%plot predictions for single growth conditions
figure
hold on
for i = 1:4
    plot(time_range, all_acc(i, :), 'o-', 'linewidth',2, 'markers', 15)
    [val, i_in] = max(all_acc(i, :));
    scatter(time_range(i_in), val, 200, 'filled', 'r', 'HandleVisibility','off')
end
ylim([60 100])
xlim([40 100])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)

%% Pairs of growth conditions
all_acc = zeros(6, length(time_range));
for i = 5:10
    for j = 1:length(time_range)
        data = [data(:, 1), data_all{i, j}];
        [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);
        all_acc(i-4, j) = max_acc;
    end
end

%plot predictions for pairs of growth conditions
figure
hold on
for i = 1:6
    plot(time_range, all_acc(i, :), 'o-', 'linewidth',2, 'markers', 15)
    [val, i_in] = max(all_acc(i, :));
    scatter(time_range(i_in), val, 200, 'filled', 'r', 'HandleVisibility','off')
end
ylim([60 100])
xlim([40 100])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)

%% Triple growth conditions
all_acc = zeros(4, length(time_range));
for i = 11:14
    for j = 1:length(time_range)
        data = [data(:, 1), data_all{i, j}];
        [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);
        all_acc(i-10, j) = max_acc;
    end
end

%plot predictions for triple growth conditions
figure
hold on
for i = 1:4
    plot(time_range, all_acc(i, :), 'o-', 'linewidth',2, 'markers', 15)
    [val, i_in] = max(all_acc(i, :));
    scatter(time_range(i_in), val, 200, 'filled', 'r', 'HandleVisibility','off')
end
ylim([60 100])
xlim([40 100])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)

%% Four growth conditions
all_acc = zeros(1, length(time_range));
for i = 15
    for j = 1:length(time_range)
        data = [data(:, 1), data_all{i, j}];
        [per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt(data, repl, fold, holdout, topk, kernelOpt);
        all_acc(1, j) = max_acc;
    end
end

%plot predictions for four growth conditions
figure
hold on
for i = 1
    plot(time_range, all_acc(i, :), 'o-', 'linewidth',2, 'markers', 15)
    [val, i_in] = max(all_acc(i, :));
    scatter(time_range(i_in), val, 200, 'filled', 'r', 'HandleVisibility','off')
end
ylim([60 100])
xlim([40 100])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)
