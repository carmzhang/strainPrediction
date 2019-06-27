%% 
% The script is for the prediction of antibiotic resistance for the library
% of clinical isolates. Here, we do 4 fold cross validation where the 
% full dataset consists of 203 strains has 4 replicates. Below are examples
% of how to run the predictions. To choose which growth condition to run 
% the prediction for, go to line 234 and pick an option described in the
% comments. To choose which antibiotic to run the prediction on, go to 
% line 194 and pick the option described in the comments. Supplementary
% Tables 1.5-1.8. This script uses the growth rate as the features.
% 
%%

close all;
clear;

close all;
clear;

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

% % anderson isolates - LB + lambda phage (MOI = 1) - 10,000x
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
datFin2 = datFin2 - min(min(datFin1))+0.0000001;
datFin3 = datFin3 - min(min(datFin1))+0.0000001;
datFin4 = datFin4 - min(min(datFin1))+0.0000001;


%alternative: convert raw data to log derivative
time = 10/60;
datFin1 = gradient(log(datFin1), time);
datFin2 = gradient(log(datFin2), time);
datFin3 = gradient(log(datFin3), time);
datFin4 = gradient(log(datFin4), time);
s = size(datFin1);

%% using only sequence isolates - comment out to use all isolates
load('resistance_WGSAnalysis_20180501.mat')
sequencedOnly = cellfun(@(x) x(4:end), sNames_WGS(1:185), 'UniformOutput', false);
sequencedNames = cellfun(@str2num, sequencedOnly);
notSeq = setdiff(myLabel, sequencedNames);
index = find(ismember(myLabel, notSeq)); %index of not sequenced
index4 = [];
for i = 1:1:size(index)
    index4 = [index4, (index(i)-1)*4+1:1:(index(i)-1)*4+4];
end
datFin1(index4, :) = [];
datFin2(index4, :) = [];
datFin3(index4, :) = [];
datFin4(index4, :) = [];
s = size(datFin1);

%load resistance information based on clinical microbiolgoy lab - true label
load('resistancesLabel.mat')

%% CHOOSE antibiotic 
% OPTIONS: 
%          1 = SAM 
%          2 = GM 
%          3 = SXT 
%          4 = CIP
%Ex: SAM is chosen
%replace following line with label = resistancesFowlerAnderson(:,2) to do
%prediction for GM
label = resistancesFowlerAnderson(:,4);%pick which antibiotic to model

label(index, :) = []; % sequenced isolates only

%% CHOOSE growth condition to do prediction on
% OPTIONS - 1: 10,000x
%           2: phage
%           3: 100x
%           4: carbenicillin
%           5: 10,000x and phage
%           6: 10,000x and 100x
%           7: 10,000x and carbenicillin
%           8: phage and 100x
%           9: phage and carbenicillin
%           10: 100x and carbenicillin
%           11: 10,000x and phage and 100x
%           12: 10,000x and phage and carbenicillin
%           13: 10,000x and 100x and carbenicillin
%           14: phage and 100x and carbenicillin
%           15: 10,000x and phage and 100x and carbenicillin
data_all = {datFin1; 
            datFin2; 
            datFin3; 
            datFin4;
            [datFin1, datFin2];
            [datFin1, datFin3];
            [datFin1, datFin4];
            [datFin2, datFin3];
            [datFin2, datFin4];
            [datFin3, datFin4];
            [datFin1, datFin2, datFin3];
            [datFin1, datFin2, datFin4];
            [datFin1, datFin3, datFin4];
            [datFin2, datFin3, datFin4];
            [datFin1, datFin2, datFin3, datFin4]
            };
%% MODIFY THIS LINE TO CHANGE GROWTH CONDITION
%Ex: prediction using 10,000x growth condition
% change the following line to change growth condition of interest
% ex - data = data_all{2} for phage condition
data = data_all{15};

%output is the optimized hyperparameters and corresponding accuracy, true
%positive rate, and true negative rate respectively
[per_kernel, idx] = clinResSVMOpt_mod(data, label);
per_kernel{1, 2}{idx, :}

disp('SAM resistance prediction using 10,000x growth condition: ')
disp('Optimized SVM hyperparameters: ')
disp(per_kernel{1, 2}{idx, 1})
disp('Accuracy: ')
disp(per_kernel{1, 2}{idx, 2})
disp('True positive rate (sensitivity): ')
disp(per_kernel{1, 2}{idx, 3})
disp('True negative rate (specificity): ')
disp(per_kernel{1, 2}{idx, 4})
