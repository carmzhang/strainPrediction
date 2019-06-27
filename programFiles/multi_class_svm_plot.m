function [predict_label,result]=multi_class_svm_plot(labeltest, labeltrain, test_feature, train_feature, N_test, N_train, kernel, cost, kkt, sigma)
%FUNCTION: returns the results of multi class svm and plots the 
%          svm boundary
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
%   predict_label - predicted label for the test set
%   result - accuracy for the test set


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
        figure
        svmOutput = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'rbf_sigma', sigma, 'kktviolationlevel', kkt, 'showplot', true);
        
    else
        svmOutput = svmtrain(train_feature,train_label, 'kernel_function', char(kernel),'boxconstraint', cost, 'kktviolationlevel', kkt, 'showplot', true);
    end
    SVMStruct{i} = svmOutput;
    
    hold on;
    for j = 1:12
        if(i == 1 && j <= 3)%if(j <=3)
            col = [255 204 51]/255;
        elseif(i == 2 && j>3 && j<=6)%elseif(j>3 && j<=6)
            col = [146 209 104]/255;
        elseif(i == 3 && j>6 && j<=9)%elseif(j>6 && j<=9)
            col = [119 238 230]/255;
        elseif(i == 4 && j>9)
            col = [245 122 202]/255;
        else
            col = [0.7 0.7 0.7];%grey for incorrect points
        end
        scatter(train_feature(j, 1), train_feature(j, 2), 600, 'MarkerFaceColor', col, 'MarkerEdgeColor', [1 1 1])
    end
    scatter(test_feature(1, 1), test_feature(1, 2), 600, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 0 0], 'LineWidth', 2)
    set(gca, 'fontsize', 26, 'FontWeight','bold','linewidth',4, 'xticklabel',[], 'yticklabel', [])
    legend('off')
    ax = gca;
    ylim([-0.1 0.15])
    xlim([-0.25 0.25])
    %ax.XColor = [186 0 10]/255;%'r'
    %ax.YColor = [186 0 10]/255;
    ax.XColor = [0 0 0]/255;
    ax.YColor = [0 0 0]/255;
    hold off;
    
    [~, F(i,:)]=svmdecision(test_feature,svmOutput);
end
[~, predict_label.svm]=max(F);
accuracy_matlab_svm_mutli=1-nnz(predict_label.svm-labeltest')/N_test;
result.svm=accuracy_matlab_svm_mutli*100;
result.distanceMatrix = F;

%temp = fitSVMPosterior(SVMStruct{1}, test_feature, labeltest)
disp(['Test accuracy, SVM, multi-class: '  num2str(result.svm),'%'])







