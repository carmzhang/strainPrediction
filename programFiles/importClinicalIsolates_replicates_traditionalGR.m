
%FUNCTION: returns the dataset for clinical isolates with all unique strains with the features
%chosen as the input parameter
% INPUT:
%   outputType - what features of the data will be in the output
%                options: 'gr', 'gc', 'AUC', 'u', 'u+AUC'
%           'gr': growth rates - derivative of the log of growth curves
%           'gc': smoothed growth curves
%           'AUC': area under the curve for the smoothed growth curves
%           'u': maximum growth rate of growth rates curve
%           'u+AUC': two features, u and AUC 
% OUTPUT: data - dataset of all unique strains with 3 biological replicates each
%                with 4 technical replicate per isolate
%                this dataset is missing strains 60-71 from the anderson library

outputType = 'gr';
seq = 'snp';%options are either snp or mlst

%import files
files = 'no antibiotic data_clean_all.xlsx';

allDat = [];
dat = xlsread(files);
%remove negative values from blanked data
dat(dat < 0) = 0;
allDat = [allDat; dat];

str = allDat(:, 1);

datFin = allDat(:, 4:end);

s = size(datFin);
%1D median filter (window = 5)
for i = 1:1:s(1)
    datFin(i, :) = conv2(datFin(i, :), ones(1, 5)/5, 'same');
end
datFin = datFin - min(min(datFin))+0.0000001;

if(ismember(outputType, {'AUC', 'u+AUC'}))
    % calculate AUC 
    AUC = trapz(datFin, 2);
    
    if(ismember(outputType, {'AUC'}))
        datFin = AUC;
    end
end

if(ismember(outputType, {'gr', 'u', 'u+AUC'}))
    %alternative: convert raw data to derivative
    time = 10/60;
    datFin = gradient(log(datFin), time);

    
    if(ismember(outputType, {'u', 'u+AUC'}))
        datFin = max(datFin, [], 2);

    end
end

if(ismember(outputType, {'u+AUC'}))
    datFin = [datFin, AUC];
end

% %% visualize
% biorep = dat(2:13, 2);
% bio = unique(biorep);
% col = {'r', 'g', 'b'};
% 
% str_unique = unique(str);
% st = ((length(str_unique)-1) * 3);
% en = st+2;
% st_e = [[0:3:st]', [2:3:en]'];
% 
% figure
% for j = 1:length(st_e)
%     ind = 0;
%     for i = st_e(j, 1):st_e(j, 2)
%         ind = ind + 1;
%         subplot(2, 8, j)
%         hold on;
%         plot(1:144, datFin((i*4 + 1):(i*4+4), :), col{ind})
%         ylabel(str_unique(j));
%     end
% end

%experimentally missing strains
missingF = [2375, 2393, 3253, 3323, 4499, 5443, 5516, 5554, 5992];
missingA = [12, 23, 26, 28, 29, 60, 61, 62, 63, 64, 65, 66, 67, 68, 70, 71];

if(ismember(seq, {'snp'}))
    %use SNP approach to cluster isolates
    load('pairwiseDist_combinedTree.mat')
    distNoRef = distPairswithout28(distPairswithout28(:,1)~=-1,:);%remove ref strain
else
    %use MLST approach to cluster isolates
    distNoRef  = xlsread('Pairwise_distance20191112_MLST.xlsx');
end

%remove experimentally missing strains
for i = 1:1:length(missingF)
    distNoRef = distNoRef(distNoRef(:,1)~=missingF(i),:);
    distNoRef = distNoRef(distNoRef(:,2)~=missingF(i),:);
end
for i = 1:1:length(missingA)
    distNoRef = distNoRef(distNoRef(:,1)~=missingA(i),:);
    distNoRef = distNoRef(distNoRef(:,2)~=missingA(i),:);
end

isolates = unique([distNoRef(:, 1); distNoRef(:, 2)]);%244 sequenced isolates
dist = zeros(length(isolates), length(isolates));

%cluster identical strains - distance = 0
identical = distNoRef(find(distNoRef(:, 3) == 0), :);
strains = cell(length(isolates), 1);
found = [];
prev = identical(end, 2);
count = 1;
for i = 0:1:length(identical)-1
    temp = identical(end-i, 1:2);
    if(temp(1, 2) == prev || ~isempty(find(strains{count, 1} == temp(1, 2))) )
        strains{count, 1} = unique([strains{count, 1}, temp]);
        found = [found, temp];
    else
        count = count + 1;
        strains{count, 1} = unique([strains{count, 1}, temp]);
    end  
    prev = temp(1, 2);
end

%strFin - each cell represents a strain cluster of isolates
strFin = strains(1:find(cellfun(@isempty,strains),1)-1, 1);
for i = 1:1:find(cellfun(@isempty,strains),1)-1
    for j = i+1:1:find(cellfun(@isempty,strains),1)-1
        if(j ~= i)
            if( ~isempty(intersect(strains{i, 1}, strains{j, 1})) && ~isempty(strFin{i, 1}))%share elements
                %strFin{i, 1} = unique([strFin{i, 1}, strains{i, 1}, strains{j, 1}]);
                strFin{i, 1} = unique([strains{i, 1}, strains{j, 1}]);
                strFin{j, 1} = [];
                
            end
            
        end
    end
end

%final list of clustered isolates
clusteredStrains = strFin(~cellfun('isempty',strFin));%clustered isolates
s = size(clusteredStrains);

%pick the first strain in the cluster for sets of clusters
pickOneC = [];
for i = 1:length(clusteredStrains)
    pickOneC = [pickOneC; clusteredStrains{i ,1}(1)];
end

%clustered isolates
allIdentical = unique([identical(:, 1); identical(:, 2)]);
allUnique = setdiff(isolates, allIdentical);%unique isolates

%set of all unique strains
uniqueStrains = [allUnique; pickOneC];

%final dataset for training and testing
data_all = [];
for i = 1:length(uniqueStrains)
    index = find(str == uniqueStrains(i));
    data_all = [data_all; str(index), repmat(i, 12, 1), datFin(index, :)];
end
data = data_all(:, 2:end);

%% strain prediction

%example 1: 3 fold cross validation with holdout (supplementary table 1.12 
% column 3) using growth rates as the features
% OUTPUT of interest: optimized - optimal hyperparameters (kernel, cost, kkt, sigma)
%                     max_acc - maximum accuracy of cross validation
%                     corresponding to the optimal hyperparameter
% (with holdout ONLY) validation_acc - accuracy of validation set
repl = 12;
fold = 2;
holdout = 1;%with holdout
topk = 1;
kernelOpt = 0;%full hyperparameter search
[per_kernel, avg_per_kernel, max_acc, optimized, validation_acc, train_acc] = clinSVMOpt_bioRep(data, repl, fold, holdout, topk, kernelOpt);

%example 2: 3 fold cross validation (supplementary table 1.12 column 2)
%using growth rates as the features
repl = 12;
fold = 3;
holdout = 0;
topk = 1;
kernelOpt = 0;
[per_kernel_2, avg_per_kernel_2, max_acc_2, optimized_2, validation_acc_2, train_acc_2] = clinSVMOpt_bioRep(data, repl, fold, holdout, topk, kernelOpt);

%Ex 1 output
disp('2 fold CV with holdout - 10,000x growth condition')
disp('Optimal SVM hyperparameters: ')
disp(optimized)
disp('CV accuracy: ')
disp(max_acc)
disp('Validation set accuracy: ')
disp(validation_acc)

%Ex 2 output
disp('3 fold CV without holdout - 10,000x growth condition')
disp('Optimal SVM hyperparameters: ')
disp(optimized_2)
disp('CV accuracy: ')
disp(max_acc_2)


