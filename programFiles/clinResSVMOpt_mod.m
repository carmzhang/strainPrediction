function [per_kernel, idx] = clinResSVMOpt_mod(data, label)
%FUNCTION: Optimizes hyperparameters of svm and returns one
%          of the optimimal parameter sets for the prediction of antibiotic
%          resistance. Uses a modified leave one out procedure.
%
% INPUT:
%   data - first column is label, other columns are features
%   label - the true resistance label
% OUTPUT:
%   per_kernel - cell with the accuracy, true positive rate, and true 
%                negative rate and predicted labels of the test
%                set for all hyperparameter sets
%  idx - index in per_kernel of the optimal parameter set

%% optimize parameter selection for svm
kernelOpt = [string('rbf')];
num_kernel = 1;
sigmaGradient = [1, 10, 20];%rbf kernel only
costGradient = [10, 100, 1000];
kktviolationlevelGradient = [0, 0.1];

% %sub-param
% kernelOpt = [string('rbf')];
% num_kernel = 1;
% sigmaGradient = [1, 20];%rbf kernel only
% costGradient = [10, 100, 1000];
% kktviolationlevelGradient = [0];

% %single param
% kernelOpt = [string('rbf')];
% num_kernel = 1;
% sigmaGradient = [20];%rbf kernel only
% costGradient = [100];
% kktviolationlevelGradient = [0.2];


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
    per_kernel{k, 2} = cell(p_size, 4);
    
    for j = 1:1:p_size
        
        s = size(data);
        m = [];
        saveScores = [];
        
        curr_param = paramAll(:,j);
        cost = curr_param(1);
        kkt = curr_param(2);
        
        %store parameter combination
        if(kernel ~= 'rbf')
            param = strcat('kernel: ', kernel, ', cost: ', num2str(cost), ', kkt: ', num2str(kkt));
        else
            sigma = curr_param(3);
            param = strcat('kernel: ', kernel, ', cost: ', num2str(cost), ', kkt: ', num2str(kkt), ', sigma: ', num2str(sigma));
        end
        
        per_kernel{k, 2}{j, 1} = param;
        
        correct = zeros(s(1)/4, 1);
        result = zeros(s(1)/4, 4);
        
        for i = 0:1:(s(1)/4) - 1
            temp = 1:1:s(1);
            temp(i*4+1:i*4+4) = [];
            train = data(temp, :);
            test = data(i*4+1:i*4+4, :);
            
            n=4 ; x=label; 
            r=repmat(x,1,n)';
            r=r(:);
            if(kernel == 'rbf')
                SVMStruct = svmtrain(train,r(temp), 'kernel_function', char(kernel),'boxconstraint', cost, 'rbf_sigma', sigma, 'kktviolationlevel', kkt, 'options', statset('MaxIter', 10000000));
            else
                SVMStruct = svmtrain(train,r(temp), 'kernel_function', char(kernel),'boxconstraint', cost, 'kktviolationlevel', kkt, 'options', optimset('MaxIter', 30000));
            end
                        F = svmclassify(SVMStruct, test);
            correct(i+1, :) = sum(F==r(i*4+1:i*4+4));
            result(i+1, :) = F';
            
            if(mod(i, 100) == 0)
                display(i)
            end
            
        end
        
        tpr = sum(result+label == 2)/sum(label == 1);
        tnr = sum(result+label == 0)/sum(label == 0);
        
        acc = sum(correct)/s(1);
        tpr = mean(tpr);%mean across 4 replicates
        tnr = mean(tnr);%mean across 4 replicates
        
        per_kernel{k, 2}{j, 2} = acc;
        per_kernel{k, 2}{j, 3} = tpr;
        per_kernel{k, 2}{j, 4} = tnr;
        
        disp(acc)
        disp(j)
        
    end %end params
    
    disp('... ...')
    disp(k)
    
end%end kernel options

all = [];
s = size(per_kernel{1, 2});
for j = 1:1:s(1)
    all = [all; mean([per_kernel{1, 2}{j, 2}, per_kernel{1, 2}{j, 3}, per_kernel{1, 2}{j, 4}])];
end

[val, idx] = max(all);
%output = per_kernel{1, 2}{idx, :};


end

