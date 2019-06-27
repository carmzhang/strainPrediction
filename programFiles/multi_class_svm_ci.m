function [out,predict_label,result, F, out2, predict_label2, result2]=multi_class_svm_ci(labeltest, labeltrain, test_feature, train_feature, N_test, N_train, kernel, cost, kkt, sigma)
%FUNCTION: returns the results of multi class svm
% INPUT:
%   labeltest - labels corresponding to the test set ( N_test by 1 matrix)
%   labeltrain - labels corrsponding to training set (N_train by 1 matrix)
%   test_feature - the test dataset (M by N_test matrix)
%   train_feature - the training dataset ( M by N_train matrix)
%   N_test - the number of samples in the test set
%   N_train -  the number of samples in the training set
%   kernel - what kernel to use for svm (ex options: 'linear', 'rbf', 'quadratic')
%   cost - svm cost hyperparameter
%   kkt - svm kkt hyperparameter
%   sigma - only a parameter for rbf kernels, if kernel is not rbf, just
%   put 0
%
% OUTPUT: 
%   out - maximum margin for each class for each sample in the test set
%   predict_label - predicted label for the test set
%   result - accuracy for the test set
%   F - svm decisions for the test set
%   out2 - maximum margin for each class for each sample in the training set
%   predict_label2 - predicted label for the training set
%   result2 - accuracy for the training set



%% 
K=max(labeltrain);

%% multi-class matlab svm
SVMStruct=cell(K);
F=zeros(K,N_test);
% train weight
for i=1:K
    train_label=ones(1,N_train);
    train_label(labeltrain==i)=-1;
    
    if(strcmp(kernel, 'rbf') == 1)
        SVMStruct{i} = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'rbf_sigma', sigma, 'kktviolationlevel', kkt);
    else
        SVMStruct{i} = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'kktviolationlevel', kkt);
    end
    [~, F(i,:)]=svmdecision(test_feature,SVMStruct{i});
    [~, F2(i, :)] = svmdecision(train_feature, SVMStruct{i});
end
out = max(F);%rows represent strain (classifiers), columns represent test set
[~, predict_label.svm]= max(F);
accuracy_matlab_svm_mutli=1-nnz(predict_label.svm-labeltest')/N_test;
result.svm=accuracy_matlab_svm_mutli*100;
result.distanceMatrix = F;

out2 = max(F2);
[~, predict_label2.svm]= max(F2);
accuracy_matlab_svm_mutli=1-nnz(predict_label2.svm-labeltrain')/N_train;
result2.svm=accuracy_matlab_svm_mutli*100;
result2.distanceMatrix = F2;

%temp = fitSVMPosterior(SVMStruct{1}, test_feature, labeltest)
disp(['Test accuracy, SVM, multi-class: '  num2str(result.svm),'%'])
disp(['Train accuracy, SVM, multi-class: '  num2str(result2.svm),'%'])







