function [result]=multi_class_svm_cv_topk(labeltest, labeltrain, test_feature, train_feature, N_test, N_train, kernel, cost, kkt, sigma, k)
% train_feature is an M by N_train matrix where N_train is the number of train data
% test_feature is an M by N_test matrix where N_test is the number of test data
% labeltrain is an N_train by 1 matrix 
% labeltest is an N_test by 1 matrix
% The 'train weight' part is slow. 


%% 
K=max(labeltrain);

%% multi-class matlab svm
SVMStruct=cell(K);
F=zeros(K,N_test);
% train weight
for i=1:K
    train_label=ones(1,N_train);
    train_label(labeltrain==i)=-1;
    
    if(kernel == 'rbf')
        SVMStruct{i} = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'rbf_sigma', sigma, 'kktviolationlevel', kkt, 'options', statset('MaxIter', 1000000));
    else
        SVMStruct{i} = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'kktviolationlevel', kkt, 'options', statset('MaxIter', 1000000));
    end
    [~, F(i,:)]=svmdecision(test_feature,SVMStruct{i});
end
[B, I] = sort(F, 'descend');
I_topk = I(1:k, :);

within_topk = 0;
for i = 1:1:N_test
    if(ismembertol(labeltest(i), I_topk(:, i)))
        within_topk = within_topk + 1;
    end
end
accuracy_matlab_svm_mutli = within_topk/N_test;

result.svm=accuracy_matlab_svm_mutli*100;
result.distanceMatrix = F;
disp(['Test accuracy, SVM, multi-class: '  num2str(result.svm),'%'])







