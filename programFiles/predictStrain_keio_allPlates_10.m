%%
% The script predicts the genetic identity of all strains in the keio
% collection. We use a 3-fold cross validation procedure to optimize the svm
% hyperparameters. All predictions are described in Supplementary Table
% 1.10. Choose what k to use in line 116.
%
% The growth dynamics data comes from this paper:
% Genome-Wide Assessment of Outer Membrane Vesicle Production in Escherichia coli
%
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

filesN = 1:2:95;
files = cell(1, length(filesN));

count = 1;
for i = 1:2:95
    files{1, count} = ['Plate', ' ', num2str(i), ' Trial 1.csv'];
    count = count + 1;
    files{1, count} = ['Plate', ' ', num2str(i), ' Trial 2.csv'];
    count = count + 1;
    files{1, count} = ['Plate', ' ', num2str(i), ' Trial 3.csv'];
    count = count + 1;
end

s = size(files);
keepStrains = [];
rawDat = [];
for i = 1:1:length(files)%s(2)
    
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
    
    %remove WT
    ind2 = strfind(names, 'WT');
    tf2 = cellfun('isempty',ind2); % true for empty cells
    ind2(tf2) = {0};
    ind2 = cell2mat(ind2);
    
    ind = or(ind, ind2);
    
    indR = find(ind > 0);
    indK = find(ind == 0);
    keepStrains = [keepStrains, names(indK)];
    
    data = gc;
    data(indR, :) = []; %remove NA
    if(mod(i, 3) == 1)
        temp{1, 1} = data;
    elseif(mod(i, 3) == 2)
        temp{2, 1} = data;
    else
        temp{3, 1} = data;
        merge = [];
        for k = 1:1:length(indK) % remove NA
        %for k = 1:1:length(data)%all data
            merge = [merge; temp{1, 1}(k, :); temp{2, 1}(k, :); ...
                temp{3, 1}(k, :)];
        end
        rawDat = [rawDat; merge];
    end
    
end

%numStrains = length(unique(keepStrains)); % remove NA and WT
numStrains = length(rawDat)/3;

%convert raw data to growth rates
s = size(rawDat);
datFin = [];
%1D moving average filter (window = 5)
for i = 1:1:s(1)
    datFin(i, :) = conv2(rawDat(i, :), ones(1, 5)/5, 'same');
end
datFin = datFin - min(min(datFin))+0.0000001;

%alternative: convert raw data to derivative
time = 10/60;
data = gradient(log(datFin), time);

% optimize parameter selection for svm
num_kernel = 1;
sigmaGradient = [1, 10, 20];%rbf kernel only
costGradient = [10, 100];
kktviolationlevelGradient = [0, 0.05];

%% CHOOSE TOP K PREDICTIONS TO LOOK AT
% OPTIONS: k = 10, 50, or 100
% Ex: k = 10
% Ex: To set k = 50 then replace line below with:
% topk = 50
topk = 10;

%split to test set 1 replicate of each strain and train set - other 3
%for all parameter pairings, split train set to train and validation set
%select highest accuracy (averaged) parameter sets

% prediction
per_kernel = cell(num_kernel, 2);%store accuracy and predictions
for k = 1:1:num_kernel
    kernel = 'rbf';
    %combination of parameters
    paramAll = combvec(costGradient, kktviolationlevelGradient, sigmaGradient);
    
    p_size = length(paramAll);
    
    per_kernel{k, 1} = kernel;
    per_kernel{k, 2} = cell(p_size, 3);
    
    for n = 1:1:p_size
        
        s = size(data);
        m = [];
        
        curr_param = paramAll(:,n);
        cost = curr_param(1);
        kkt = curr_param(2);
        
        %store parameter combination
        sigma = curr_param(3);
        param = strcat('kernel: ', kernel, ', cost: ', num2str(cost), ', kkt: ', num2str(kkt), ', sigma: ', num2str(sigma));
        per_kernel{k, 2}{n, 1} = param;
        
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
            [result]=multi_class_svm_cv_topk(test_label,train_label,test_feature,train_feature,N_test,N_train, kernel, cost, kkt, sigma, topk);
            % result.svm is the accuracy of svm
            
            m = [m, result.svm];%store accuracy
            disp(i)
            
        end
        
        per_kernel{k, 2}{n, 2} = m;
        
        disp('...')
        disp(n)
        save('all_plates_top10_cv_20190516.mat')
        
    end %end params
    
    disp('... ...')
    disp(k)
    
end%end kernel options

%take average across 3 folds and find optimal parameter set
avg_per_kernel = cell(num_kernel, 2);
for i = 1:1:num_kernel
    avg_per_kernel{i, 1} = per_kernel{i, 1};
    s = size(per_kernel{i, 2});
    all = [];
    for o = 1:1:s(1)
        all = [all, mean(per_kernel{i, 2}{o, 2})];
    end
    avg_per_kernel{i, 2} = all;
end

optimal_out = cell(num_kernel, 2);
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

save('all_plates_top10_cv_20190516.mat')

%1x2 cell: (1) optimized hyperparameter (2) classification accuracy
optimal_out



