%%
% Script plots figure 3c in the main text:
% (1) The histogram of the maximum margins corresponding to correctly and
%     incorrectly classified samples
% (2) The plot with the logistic regression line and 95% confidence 
%     interval to estimate the probability of a prediction being correct
%     
%%

close all;
clear;

data = importClinicalIsolates('gr');

s = size(data);
saveScores = [];
margin = [];
m = [];
saveScores2 = [];
margin2 = [];
m2 = [];
F_final = cell(4, 1);
for i = 1:1:4
    
    temp = 1:1:s(1);
    temp(i:4:end) = [];
    trainX = data(temp, 2:99); %10,000x
    testX = data(i:4:end, 2:99);

    N_test = 203;% number of testing samples
    N_train = 203*3;% number of training samples
    
    test_label = data(i:4:end, 1);
    train_label = data(temp, 1);
    test_feature = testX;
    train_feature = trainX;
    
    %%
    [f, predict_label,result, F_all, f2, predict_label2, result2]=multi_class_svm_ci(test_label,train_label,test_feature,train_feature,N_test,N_train, 'rbf', 100, 0, 10);
    
    F_final{i, 1} = F_all;%rows represent strain (classifiers), columns represent test set
    
    %for test set
    m = [m, result.svm];%store accuracy
    saveScores = [saveScores, predict_label.svm'];
    margin = [margin, f'];
    
    % for train set
    m2 = [m2, result2.svm];%store accuracy
    saveScores2 = [saveScores2, predict_label2.svm'];
    margin2 = [margin2, f2'];
    
    disp(i)
    
end

misclass = cell(4, 2);%column 1 is index, column 2 is margin
all_miscl_margin = [];
correct = cell(4, 2);
all_correct_margin = [];
for i = 1:1:4
    temp = find(saveScores(:, 1) - [test_label] ~= 0);
    misclass{i, 1} = temp;
    misclass{i, 2} = margin(temp, i);
    all_miscl_margin = [all_miscl_margin; margin(temp, i)];
    
    temp2 = find(saveScores(:, 1) - [test_label] == 0);
    correct{i, 1} = temp2;
    correct{i, 2} = margin(temp2, i);
    all_correct_margin = [all_correct_margin; margin(temp2, i)];
end

corr = repmat(1, [length(all_correct_margin), 1]);
inc = repmat(0, [length(all_miscl_margin), 1]);


bins = floor(sqrt(203*4));

[series1, centers] = hist(reshape(margin, [812, 1]), bins);

x_all = [all_correct_margin; all_miscl_margin];
y_all = [corr; inc];

dat_cl = [x_all, y_all];

[b, dev, stats] = glmfit(x_all,y_all, 'binomial', 'link', 'logit');
[yfit ,lo, hi] = glmval(b,x_all, 'logit', stats, 'size', 1, 'Confidence', 0.95, 'simultaneous' , true);
fit_out = [x_all, yfit, lo, hi];
fit_out = sortrows(fit_out, 1);

%plot probability of prediction being correct
figure
hold on
plot([0 0], [0 1], '--', 'Color', [0 0 0], 'linewidth', 2)
ciplot(fit_out(:, 2)-fit_out(:, 3), fit_out(:, 2)+fit_out(:, 4), fit_out(:, 1), [206/255 208/255 206/255])
plot(fit_out(:, 1), fit_out(:, 2),'-','Color', [69/255 72/255 81/255] ,'LineWidth',5)
plot(fit_out(:, 1), fit_out(:, 2)-fit_out(:, 3), '--', 'Color', [141/255 153/255 174/255])
plot(fit_out(:, 1), fit_out(:, 2)+fit_out(:, 4), '--', 'Color', [141/255 153/255 174/255])
xlim([min(fit_out(:, 1)), max(fit_out(:, 1))])
set(gca, 'fontsize', 48, 'FontWeight','bold','linewidth',4)

%histogram of correct and incorrect predictions per maximum margin
figure
histf(all_correct_margin, centers, 'EdgeColor', [34 112 39]/255, 'FaceColor', [222 255 219]/255, 'linewidth', 2)
box off;
hold on;
histf(all_miscl_margin, centers, 'EdgeColor', [217 30 54]/255, 'FaceColor', [255 196 196]/255, 'linewidth', 2)
set(gca, 'fontsize', 48, 'FontWeight','bold','linewidth',4)
xlim([min(fit_out(:, 1)), max(fit_out(:, 1))])
ylim([0 80])
legend('boxoff')
