function roc_all = resistance_roc(data_all, antibiotic, param, index)
%FUNCTION: plots roc curve 
% INPUT:
%   data_all - cell containing all combinations of growth conditions
%   antibiotic - which antibiotic to do prediction for
%          OPTIONS: 1 = SAM   2 = GM   3 = SXT   4 = CIP
%   param - cell containing optimized parameter sets for all combinations
%           of growth conditions
%   index - non sequenced isolates
%
% OUTPUT: roc curve for all combinations of growth conditions
%   roc_all - data used to plot roc curve for all growth conditions

%svm prediction of resistance - 1:SAM, 2:GM, 3: SXT, 4:CIP
load('resistancesLabel.mat')
label = resistancesFowlerAnderson(:, antibiotic);%pick which antibiotic to model
label(index, :) = []; %%%% comment out for all isolates

s = size(data_all{1});

roc_all = cell(length(data_all), 1);
for j = 1:length(data_all)
    %data: 1. 10,000x 2. phage 3. 100x, 4. carb
    data = data_all{j};
    
    
    %% optimized parameters for svm
    kernel = string('rbf');
    cost = param(j, 1);
    kkt = param(j, 2);
    sigma = param(j, 3);
    
    %% prediction
%     correct = zeros(s(1)/4, 1);
%     result = zeros(s(1)/4, 4);
    margin = zeros(s(1)/4, 4);
    
    for i = 0:1:(s(1)/4) - 1
        temp = 1:1:s(1);
        temp(i*4+1:i*4+4) = [];
        train = data(temp, :);
        test = data(i*4+1:i*4+4, :);
        
        n=4 ; x=label;
        r=repmat(x,1,n)';
        r=r(:);
        
        if(kernel == 'rbf')
            SVMStruct = svmtrain(train,r(temp), 'kernel_function', char(kernel),'boxconstraint', cost, 'rbf_sigma', sigma, 'kktviolationlevel', kkt, 'options', optimset('MaxIter', 30000));
        else
            SVMStruct = svmtrain(train,r(temp), 'kernel_function', char(kernel),'boxconstraint', cost, 'kktviolationlevel', kkt, 'options', optimset('MaxIter', 30000));
        end
        
%         F = svmclassify(SVMStruct, test);
%         correct(i+1, :) = sum(F==r(i*4+1:i*4+4));
%         result(i+1, :) = F';
        
        [~, marg]=svmdecision(test,SVMStruct);
        margin(i+1, :) = marg';
        
    end
%     tpr = sum(result+label == 2)/sum(label == 1);
%     tnr = sum(result+label == 0)/sum(label == 0);
%     
%     acc = sum(correct)/s(1);
%     tpr = mean(tpr);%mean across 4 replicates
%     tnr = mean(tnr);%mean across 4 replicates
    
    %negative margin value corresponds to 1 and positive margin value
    %corresponds to 0
    %FPR = 1-TNR vs TPR is ROC
    thresholds = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];
    
    result_roc = cell(length(thresholds), 1);
    for i = 1:length(thresholds)
        
        cur = thresholds(i);
        margin_temp = margin;
        idx = find(margin_temp >= cur);
        margin_temp(idx) = 11;
        idx = find(margin_temp < cur);
        margin_temp(idx) = 1;
        idx = find(margin_temp == 11);
        margin_temp(idx) = 0;
        result_roc{i, 1} = margin_temp;
    end
    
    roc = zeros(length(thresholds), 2);
    for i = 1:length(thresholds)
        tpr = sum(result_roc{i}+label == 2)/sum(label == 1);
        tnr = sum(result_roc{i}+label == 0)/sum(label == 0);
        
        tpr = mean(tpr);%mean across 4 replicates
        tnr = mean(tnr);%mean across 4 replicates
        
        roc(i, :) = [tpr, tnr];
    end
    roc_all{j, 1} = roc;
    disp(j)
end

% plot roc curve for all growth conditions - optimized svm parameters

% %4, 6, 4, 1
  
% %attempt 1 colors
% col = [linspace(255, 120, 4)', linspace(0, 0, 4)', linspace(255, 100, 4)';
%        linspace(0, 60, 6)', linspace(255, 179, 6)', linspace(0, 0, 6)';
%        linspace(6, 6, 4)', linspace(239, 20, 4)', linspace(227, 227, 4)';
%        0, 0, 0]/255;
   
%RGB - order 
col = [linspace(255, 120, 4)', linspace(114, 0, 4)', linspace(255, 100, 4)';%R
       linspace(215, 47, 6)', linspace(237, 142, 6)', linspace(151, 48, 6)';%G
       linspace(211, 127, 4)', linspace(243, 183, 4)', linspace(238, 222, 4)';%B
       0, 0, 0]/255;
   
col = [linspace(255, 120, 4)', linspace(206, 0, 4)', linspace(255, 100, 4)';%R
       linspace(215, 47, 6)', linspace(237, 142, 6)', linspace(151, 48, 6)';%G
       linspace(211, 127, 4)', linspace(243, 183, 4)', linspace(238, 222, 4)';%B
       0, 0, 0]/255;

if(size(data_all, 1) == 3)
   col = col([1, 3, 6], :);
end
%overlay on plot resistance prediction from wgs
[roc_all_resfinder, roc_all_card] = resistance_roc_WGS(antibiotic);

   
figure
hold on
for i = 1:length(data_all)
    roc = roc_all{i, 1};
    plot(1-roc(:, 2), roc(:, 1), '-o', 'linewidth', 3.5, 'color', col(i, :), 'markerfacecolor', col(i, :))
end
plot(0:0.1:1, 0:0.1:1, '--k', 'linewidth', 3)
plot(1-roc_all_resfinder(:, 2), roc_all_resfinder(:, 1), '--o', 'linewidth', 3.5, 'color', 'r', 'markerfacecolor', 'r')
plot(1-roc_all_card(:, 2), roc_all_card(:, 1), ':o', 'linewidth', 3.5, 'color', 'r', 'markerfacecolor', 'r')
%legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15')
%legend boxoff
set(gca, 'fontsize', 40, 'FontWeight','Bold', 'linewidth', 6)



end

