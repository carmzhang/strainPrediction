%%
% The script predicts the genetic identity of subsets of the keio
% collection. We use a 3-fold cross validation procedure to optimize the svm
% hyperparameters. All predictions are described in Supplementary Table
% 1.9. Choose which subset (plate) of the library to do the prediction on
% in line 22. This script uses the derivative of growth curves as the features.
%
% The growth dynamics data comes from this paper:
% Genome-Wide Assessment of Outer Membrane Vesicle Production in Escherichia coli
%     Adam J. Kulp,
%     Bo Sun,
%     Teresa Ai,
%     Andrew J. Manning,
%     Nichole Orench-Rivera,
%     Amy K. Schmid,
%     Meta J. Kuehn
%%

close all;
clear;

%% CHOOSE SUBSET TO PREDICT GENETIC IDENTITY
% OPTIONS: 1 = plate 5
%          2 = plate 7
%          3 = plate 9
%          4 = plate 11
%          5 = plate 13
%          6 = plate 15
%          7 = plate 17
%          8 = plate 19
files_all = {'Plate 5 Trial 1.csv', 'Plate 5 Trial 2.csv', 'Plate 5 Trial 3.csv';
             'Plate 7 Trial 1.csv', 'Plate 7 Trial 2.csv', 'Plate 7 Trial 3.csv';
             'Plate 9 Trial 1.csv', 'Plate 9 Trial 2.csv', 'Plate 9 Trial 3.csv';
             'Plate 11 Trial 1.csv', 'Plate 11 Trial 2.csv', 'Plate 11 Trial 3.csv';
             'Plate 13 Trial 1.csv', 'Plate 13 Trial 2.csv', 'Plate 13 Trial 3.csv';
             'Plate 15 Trial 1.csv', 'Plate 15 Trial 2.csv', 'Plate 15 Trial 3.csv';
             'Plate 17 Trial 1.csv', 'Plate 17 Trial 2.csv', 'Plate 17 Trial 3.csv';
             'Plate 19 Trial 1.csv', 'Plate 19 Trial 2.csv', 'Plate 19 Trial 3.csv';
            };
%% EX: choose plate 5
% if want to change to plate 7 replace line below with:
% files = files_all(2, :);
files = files_all(1, :);

meanAcc = [];
s = size(files);
keepStrains = [];
rawDat = [];
for i = 1:1:s(2)
    
    if(mod(i, 3) == 1)
        temp = cell(3, 1);
    end
    t5 = readtable(files{1, i});
    time = table2array(t5(:, 1));
    gc = table2array(t5(:, 2:end))';
    
    names = t5.Properties.VariableNames(2:end);
    
    %remove NA
    ind = strfind(names, 'x_N_A');
    tf = cellfun('isempty',ind); % true for empty cells
    ind(tf) = {0};
    ind = cell2mat(ind);
    indR = find(ind == 1);
    indK = find(ind == 0);
    keepStrains = [keepStrains, names(indK)];
    
    data = gc;
    data(indR, :) = []; %to remove NA
    if(mod(i, 3) == 1)
        temp{1, 1} = data;
    elseif(mod(i, 3) == 2)
        temp{2, 1} = data;
    else
        temp{3, 1} = data;
        merge = [];
        for k = 1:1:length(indK) %to remove NA
        %for k = 1:1:length(data)%all data
            merge = [merge; temp{1, 1}(k, :); temp{2, 1}(k, :); ...
                temp{3, 1}(k, :)];
        end
        rawDat = [rawDat; merge];
    end
    
end

%numStrains = length(unique(keepStrains)); % remove NA
numStrains = length(rawDat)/3;

%convert raw data to growth rates
s = size(rawDat);
datFin = [];
%1D median filter (window = 3)
for i = 1:1:s(1)
    datFin(i, :) = medfilt1(rawDat(i, :));
end

%alternative: convert raw data to derivative
newDat = [];
for i = 1:1:s(1)
    x1 = time(1:end-1);
    x2 = time(2:end);
    y1 = datFin(i, 1:end-1)';
    y2 = datFin(i, 2:end)';
    slopes = (y2 - y1) ./ (x2 - x1);
    newDat(i, :) = slopes';
end
data = newDat;


%% optimize parameter selection for svm
kernelOpt = [string('linear'), string('quadratic'), string('rbf')];
num_kernel = 3;
sigmaGradient = [1, 10, 20];%rbf kernel only
costGradient = [10, 100, 1000];
kktviolationlevelGradient = [0, 0.05, 0.1];

%split to test set 1 replicate of each strain and train set - other 3
%for all parameter pairings, split train set to train and validation set
%select highest accuracy (averaged) parameter sets

topk = 10;

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
        
        s = size(data);
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
        
        %% prediction - SVM
        data(isnan(data))= 0;%replace nan with 0
        s = size(data);
        m = [];
        saveScores = [];
        for i = 1:1:3
            
            temp = 1:1:s(1);
            temp(i:3:end) = [];
            trainX = data(temp, :);
            testX = data(i:3:end, :);
            
            N_test = numStrains;% number of testing samples
            N_train = numStrains*2;% number of training samples
            K = numStrains;%number of classes
            M = s(2);%dimension of features
            
            x=(1:N_test)';
            r=repmat(x,1,2)';
            trainY = r(:);
            
            test_label = [1:1:numStrains]';
            testY = test_label;
            train_label = trainY;
            test_feature = testX;
            train_feature = trainX;
            
            %%
            %[predict_label,result]=multi_class_svm_cv(test_label,train_label,test_feature,train_feature,N_test,N_train, kernel, cost, kkt, sigma);
            [result]=multi_class_svm_cv_topk(test_label,train_label,test_feature,train_feature,N_test,N_train, kernel, cost, kkt, sigma, topk);
            % predict_label.svm is the predicted labels based on Matlab default SVM
            % result.svm is the accuracy of svm
            
            m = [m, result.svm];%store accuracy
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

%take average across 3 folds and find optimal parameter set
avg_per_kernel = cell(3, 2);
for i = 1:1:num_kernel
    avg_per_kernel{i, 1} = per_kernel{i, 1};
    s = size(per_kernel{i, 2});
    all = [];
    for j = 1:1:s(1)
        all = [all, mean(per_kernel{i, 2}{j, 2})];
    end
    avg_per_kernel{i, 2} = all;
end

optimal_out = cell(1, 2);
cur_max = 0;
for i = 1:1:num_kernel
    [val, idx] = max(avg_per_kernel{i, 2});
    per_kernel{i, 2}{idx, 1};%optimal parameter set
    
    if(cur_max < val)
       cur_max  = val;
       optimal_out{1, 1} = per_kernel{i, 2}{idx, 1};%optimal parameter set
       optimal_out{1, 2} = val;
    end
end

%1x2 cell: (1) optimized hyperparameter (2) classification accuracy
optimal_out
