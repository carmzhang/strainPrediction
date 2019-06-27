function [saveScores] = singlePredict_mod(data)
%FUNCTION: Runs 4 fold cross validation without holdout (meant for clinical
%          isolates) and returns the predicted margin values.
% INPUT:
%   data - first column is label, other columns are features
%
% OUTPUT:
%   saveScores - margin values for all replicates when they are in the test
%   set


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
    
    N_test = size(testX, 1);%203;% number of testing samples
    N_train = size(trainX, 1);%203*3;% number of training samples
    
    test_label = data(i:4:end, 1);
    train_label = data(temp, 1);
    test_feature = testX;
    train_feature = trainX;
    
    %%
    [f, predict_label,result, F_all, f2, predict_label2, result2]=multi_class_svm_ci_mod(test_label,train_label,test_feature,train_feature,N_test,N_train, 'rbf', 100, 0.2, 10);
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


end

