function [data] = importClinicalIsolates_traditionalGR(outputType)
%FUNCTION: returns the dataset for clinical isolates with all unique strains with the features
%chosen as the input parameter
% INPUT:
%   outputType - what features of the data will be in the output
%                options: 'gr', 'gc', 'AUC', 'u', 'u+AUC'
%           'gr': growth rates - derivative of log of  growth curves 
%           'gc': smoothed growth curves
%           'AUC': area under the curve for the smoothed growth curves
%           'u': maximum value of 'gr' curve
%           'u+AUC': two features, u and AUC 
% OUTPUT: data - dataset of all unique strains with 4 replicates per strain

%% load fowler data
load('strainNames.mat')

num = 13;
files = cell(num, 1);
filesA = cell(3, 1);

% fowler isolates - LB - 10,000x 
files{1, 1} = '20170715_Fowler_1-23_LB_10000.mat';
files{2, 1} = '20170715_Fowler_24-46_LB_10000.mat';
files{3, 1} = '20170716_Fowler_47-69_LB_10000.mat';
files{4, 1} = '20170716_Fowler_70-92_LB_10000.mat';
files{5, 1} = '20170717_Fowler_93-106_LB_10000.mat';
files{6, 1} = '20170719_Fowler_107-129_LB_10000.mat';
files{7, 1} = '20170720_Fowler_130-149_LB_10000.mat';
files{8, 1} = '20170721_Fowler_150-172_LB_10000.mat';
files{9, 1} = '20170722_Fowler_173-195_LB_10000.mat';
files{10, 1} = '20170722_Fowler_196-218_LB_10000.mat';
% anderson isolates - LB - 10,000x
filesA{1, 1} = '20170805_Anderson_1-20_LB_10000.mat';
filesA{2, 1} = '20170806_Anderson_21-40_LB_10000.mat';
filesA{3, 1} = '20170806_Anderson_41-59_LB_10000.mat';

allDat = [];
for i = 1:1:10
    load(files{i, 1})
    %remove negative values from blanked data
    datFin(datFin < 0) = 0;
    allDat = [allDat; datFin];
end 
r = repmat(myLabelLB,1,4)'; %order of isolates
lab = r(:);
%reorder data to match others
allDat = [lab, allDat];
datFin1 = sortrows(allDat);
datFin1 = datFin1(:, 2:end);

allDat = [];
for i = 1:1:3
    load(filesA{i, 1})
    %remove negative values from blanked data
    datFin(datFin < 0) = 0;
    allDat = [allDat; datFin];
end 
datFin1 = [datFin1; allDat];

%fowler isolates - LB 10,000x + lambda phage MOI = 1
files{1, 1} = '20170807_Fowler_1-23_LB_lambda_MOI-1_10000.mat';
files{2, 1} = '20170807_Fowler_24-46_LB_lambda_MOI-1_10000.mat';
files{3, 1} = '20170808_Fowler_47-69_LB_lambda_MOI-1_10000.mat';
files{4, 1} = '20170808_Fowler_70-92_LB_lambda_MOI-1_10000.mat';
files{5, 1} = '20170809_Fowler_93-115_LB_lambda_MOI-1_10000.mat';
files{6, 1} = '20170809_Fowler_116-138_LB_lambda_MOI-1_10000.mat';
files{7, 1} = '20170809_Fowler_139-161_LB_lambda_MOI-1_10000.mat';
files{8, 1} = '20170810_Fowler_162-184_LB_lambda_MOI-1_10000.mat';
files{9, 1} = '20170810_Fowler_185-207_LB_lambda_MOI-1_10000.mat';
files{10, 1} = '20170810_Fowler_208-218_LB_lambda_MOI-1_10000.mat';
% anderson isolates - LB + lambda phage (MOI = 1) - 10,000x
files{11, 1} = '20170811_Anderson_1-23_LB_lambda_MOI-1_10000.mat';
files{12, 1} = '20170811_Anderson_24-46_LB_lambda_MOI-1_10000.mat';
files{13, 1} = '20170811_Anderson_47-59_LB_lambda_MOI-1_10000.mat';

allDat = [];
for i = 1:1:num
    load(files{i, 1})
    %remove negative values from blanked data
    datFin(datFin < 0) = 0;
    allDat = [allDat; datFin];
end 
datFin2 = allDat;

%fowler isolates - LB - 100x
files{1, 1} = '20170818_Fowler_1-23_LB_100.mat';
files{2, 1} = '20170818_Fowler_24-46_LB_100.mat';
files{3, 1} = '20170818_Fowler_47-69_LB_100.mat';
files{4, 1} = '20170819_Fowler_70-92_LB_100.mat';
files{5, 1} = '20170819_Fowler_93-115_LB_100.mat';
files{6, 1} = '20170819_Fowler_116-138_LB_100.mat';
files{7, 1} = '20170820_Fowler_139-161_LB_100.mat';
files{8, 1} = '20170820_Fowler_162-184_LB_100.mat';
files{9, 1} = '20170820_Fowler_185-207_LB_100.mat';
files{10, 1} = '20170821_Fowler_208-218_LB_100.mat';
%anderson isolates - LB - 100x
files{11, 1} = '20170821_Anderson_1-23_LB_100.mat';
files{12, 1} = '20170821_Anderson_24-46_LB_100.mat';
files{13, 1} = '20170821_Anderson_47-59_LB_100.mat';

allDat = [];
for i = 1:1:num
    load(files{i, 1})
    %remove negative values from blanked data
    datFin(datFin < 0) = 0;
    allDat = [allDat; datFin];
end 
datFin3 = allDat;

%fowler isolates - LB - carbenicillin (5 ug/ml) - 10,000x
files{1, 1} = '20170903_Fowler_1-23_LB_Carb-5ugml_10000.mat';
files{2, 1} = '20170904_Fowler_24-46_LB_Carb-5ugml_10000.mat';
files{3, 1} = '20170904_Fowler_47-69_LB_Carb-5ugml_10000.mat';
files{4, 1} = '20170905_Fowler_70-92_LB_Carb-5ugml_10000.mat';
files{5, 1} = '20170906_Fowler_93-115_LB_Carb-5ugml_10000.mat';
files{6, 1} = '20170909_Fowler_116-138_LB_Carb-5ugml_10000.mat';
files{7, 1} = '20170907_Fowler_139-161_LB_Carb-5ugml_10000.mat';
files{8, 1} = '20170908_Fowler_162-184_LB_Carb-5ugml_10000.mat';
files{9, 1} = '20170908_Fowler_185-207_LB_Carb-5ugml_10000.mat';
files{10, 1} = '20170909_Fowler_208-218_LB_Carb-5ugml_10000.mat';
%anderson isolates - LB - carbenicillin (5 ug/ml) - 10,000x
files{11, 1} = '20170910_Anderson_1-23_LB_Carb-5ugml_10000.mat';
files{12, 1} = '20170910_Anderson_24-46_LB_Carb-5ugml_10000.mat';
files{13, 1} = '20170911_Anderson_47-59_LB_Carb-5ugml_10000.mat';

allDat = [];
for i = 1:1:num
    load(files{i, 1})
    %remove negative values from blanked data
    datFin(datFin < 0) = 0;
    allDat = [allDat; datFin];
end 
datFin4 = allDat;

s = size(datFin1);
%1D moving average filter (window = 5)
for i = 1:1:s(1)
    datFin1(i, :) = conv2(datFin1(i, :), ones(1, 5)/5, 'same');
    datFin2(i, :) = conv2(datFin2(i, :), ones(1, 5)/5, 'same');
    datFin3(i, :) = conv2(datFin3(i, :), ones(1, 5)/5, 'same');
    datFin4(i, :) = conv2(datFin4(i, :), ones(1, 5)/5, 'same');
end

datFin1 = datFin1 - min(min(datFin1))+0.0000001;
datFin2 = datFin2 - min(min(datFin2))+0.0000001;
datFin3 = datFin3 - min(min(datFin3))+0.0000001;
datFin4 = datFin4 - min(min(datFin4))+0.0000001;

% s = size(datFin1);
% %1D median filter (window = 3)
% for i = 1:1:s(1)
%     datFin1(i, :) = medfilt1(datFin1(i, :));
%     datFin2(i, :) = medfilt1(datFin2(i, :));
%     datFin3(i, :) = medfilt1(datFin3(i, :));
%     datFin4(i, :) = medfilt1(datFin4(i, :));
% end
% datFin1 = datFin1 - min(min(datFin1))+0.0000001;
% datFin2 = datFin2 - min(min(datFin2))+0.0000001;
% datFin3 = datFin3 - min(min(datFin3))+0.0000001;
% datFin4 = datFin4 - min(min(datFin4))+0.0000001;


if(ismember(outputType, {'AUC', 'u+AUC'}))
    % calculate AUC 
    AUC1 = trapz(datFin1, 2);
    AUC2 = trapz(datFin2, 2);
    AUC3 = trapz(datFin3, 2);
    AUC4 = trapz(datFin4, 2);
    
    if(ismember(outputType, {'AUC'}))
        datFin1 = AUC1;
        datFin2 = AUC2;
        datFin3 = AUC3;
        datFin4 = AUC4;
    end
end

if(ismember(outputType, {'gr', 'u', 'u+AUC'}))
    %alternative: convert raw data to derivative
    time = 10/60;
    datFin1 = gradient(log(datFin1), time);
    datFin2 = gradient(log(datFin2), time);
    datFin3 = gradient(log(datFin3), time);
    datFin4 = gradient(log(datFin4), time);
    
    if(ismember(outputType, {'u', 'u+AUC'}))
        datFin1 = max(datFin1, [], 2);
        datFin2 = max(datFin2, [], 2);
        datFin3 = max(datFin3, [], 2);
        datFin4 = max(datFin4, [], 2);
    end
end

if(ismember(outputType, {'u+AUC'}))
    datFin1 = [datFin1, AUC1];
    datFin2 = [datFin2, AUC2];
    datFin3 = [datFin3, AUC3];
    datFin4 = [datFin4, AUC4];
end

datFin = [datFin1, datFin2, datFin3, datFin4];

%exp strains for fowler and anderson that are missing
missingF = [2375, 2393, 3253, 3323, 4499, 5443, 5516, 5554, 5992];
missingA = [12, 23, 26, 28, 29];

%load tree
load('pairwiseDist_combinedTree.mat')

%remove reference strain
distNoRef = distPairswithout28(distPairswithout28(:,1)~=-1,:);
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

%clustered isolates
allIdentical = unique([identical(:, 1); identical(:, 2)]);
allUnique = setdiff(isolates, allIdentical);%unique isolates

%reorder data and pull out 4 random replicates for clustered strains
%myLabel is order of isolates for experiments in fowler library
%myLabelA is order of isolates for experiments in anderson library (isolate ID from isolate 2 list)
myLabelA = [1:7,9, 10, 13, 16:21, 24, 25, 27, 30:55, 57:68, 70:71]';

%all clusters
s = size(clusteredStrains);
allStrains = cell(length(allUnique) + s(1), 1);%clustered strains - sequencing
s2 = size(allStrains);
j = 1;
for i = 1:1:s2(1)
    if(i <= s(1))
        allStrains{i, 1} = clusteredStrains{i, 1};
    else
        allStrains{i, 1} = [allUnique(j)];
        j = j + 1;
    end
end

allLabels = [myLabel; myLabelA];%experimental isolates
reorderDat = cell(s2(1), 1);
count = 0;
notSeq = [];
for i = 0:1:length(allLabels)-1
    
    cur = allLabels(i+1);
    c = cellfun(@(x)(ismember(cur,x)), allStrains,'UniformOutput',false);
    [YS,~] = find(reshape([c{:}],numel(cur),[])');
    
    if(~isempty(YS))
        reorderDat{YS, 1} = [reorderDat{YS, 1}; datFin(i*4+1:i*4+4, :)];
    else
        count = count + 1;
        notSeq = [notSeq; cur];
    end
end

%% data for training + testing
finalDat = reorderDat(~cellfun('isempty',reorderDat));%clustered isolates tc
s = size(finalDat);
data = [];
for i = 1:1:s(1)
    cur = finalDat{i, 1};
    sub = size(cur);
    if(sub(1) == 4)
       data = [data; i, cur(1, :); i, cur(2, :); i, cur(3, :) ...
           ; i, cur(4, :)];
    else
        temp = randperm(sub(1));
        newOrd = cur(temp(1:4), :);
        data = [data; i, newOrd(1, :); i, newOrd(2, :); ...
            i, newOrd(3, :); i, newOrd(4, :)];
    end
end

end

