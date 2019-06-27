%% 
% The script plots the left subplot in Figure 4b and all subplots for
% Figure 4c using the growth rate as features.
% (1) Figure 4b - correlation between genetic and phenotypic distance
% (2) Figure 4c -  correlation between phenotypic distances (10,000x vs 
%     phage) 
% (3) Figure 4c -  correlation between phenotypic distances (10,000x vs 
%     carbenicillin) 
% (4) Figure 4c -  correlation between phenotypic distances (phage vs 
%     carbenicillin) 
%
%% phenotype vs genotype mapping

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
datFin2 = datFin2 - min(min(datFin2))+0.0000001;
datFin3 = datFin3 - min(min(datFin3))+0.0000001;
datFin4 = datFin4 - min(min(datFin4))+0.0000001;
%convert raw data to derivative
time = 10/60;
datFin1 = gradient(log(datFin1), time);
datFin2 = gradient(log(datFin2), time);
datFin3 = gradient(log(datFin3), time);
datFin4 = gradient(log(datFin4), time);
s = size(datFin1);

%label order
myLabelA = [1:7,9, 10, 13, 16:21, 24, 25, 27, 30:55, 57:68, 70:71]';
myLabel;
label_all = [myLabel; myLabelA];

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

r=repmat(label_all,1,4)';
label_fin = r(:);
data = [label_fin, datFin1, datFin2, datFin3, datFin4];

data_m = reshape(mean(reshape(data,4,[])),[],397);%mean of 4 replicates

%resistances
load('resistancesLabel.mat')
data_r = [label_all, resistancesFowlerAnderson];

pair_wgs = zeros(length(isolates), length(isolates));

one = 2:100;
two = 101:199;
three = 200:298;
four = 299:397;

pair_phenotype_m = zeros(length(isolates), length(isolates));
pair_phenotype_1m = zeros(length(isolates), length(isolates));
pair_phenotype_2m = zeros(length(isolates), length(isolates));
pair_phenotype_3m = zeros(length(isolates), length(isolates));
pair_phenotype_4m = zeros(length(isolates), length(isolates));

for i = 1:length(isolates)
    
    i1 = isolates(i);
    
    for j = 1:length(isolates)
        i2 = isolates(j);
        if(i1 ~= i2)
            idx1_m = find(data_m(:, 1) == i1);
            idx2_m = find(data_m(:, 1) == i2);
            pair_phenotype_m(i, j) = pdist([data_m(idx1_m, 2:end); data_m(idx2_m, 2:end)]);
            
            pair_phenotype_1m(i, j) = pdist([data_m(idx1_m, one); data_m(idx2_m, one)]);
            pair_phenotype_2m(i, j) = pdist([data_m(idx1_m, two); data_m(idx2_m, two)]);
            pair_phenotype_3m(i, j) = pdist([data_m(idx1_m, three); data_m(idx2_m, three)]);
            pair_phenotype_4m(i, j) = pdist([data_m(idx1_m, four); data_m(idx2_m, four)]);
            
            if(~isempty(intersect(find(distNoRef(:, 1) == i1),find(distNoRef(:, 2) == i2))))
                row_in = distNoRef(intersect(find(distNoRef(:, 1) == i1),find(distNoRef(:, 2) == i2)), 3);
                pair_wgs(i, j) = row_in;
            else
                row_in = distNoRef(intersect(find(distNoRef(:, 2) == i1),find(distNoRef(:, 1) == i2)), 3);
                pair_wgs(i, j) = row_in;
            end
        end
        
    end
end

result_noidentical = [];
result_noidentical_m = [];
for i = 1:length(isolates)
    for j = 1:length(isolates)
        if(i ~= j)
            result_noidentical_m = [result_noidentical_m; pair_wgs(i, j), pair_phenotype_m(i, j), pair_phenotype_1m(i, j), pair_phenotype_2m(i, j), pair_phenotype_3m(i, j), pair_phenotype_4m(i, j)];%3-6 is gc condition 1-4
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 - density plot
% -  mean of 4 replicates - 1.(3) 10,000x, 2. (4) phage 3. (5) 100x 4. (6)
% carbenicillin

nbins = 8;
%wgs vs 10,000x
densityplot(result_noidentical_m(:, 1), result_noidentical_m(:,3), [50, 50])
hold on;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
[N, edges, bin] = histcounts(result_noidentical_m(:, 1), nbins);
x = [];
y = [];
sd = [];
for i = 1:1:nbins
    idx = find(bin == i);
    y = [y mean(result_noidentical_m(idx, 3))];
    x = [x mean([edges(i) edges(i+1)])];
    sd = [sd std(result_noidentical_m(idx, 3))];
end
errorbar(x, y, sd, 'linewidth', 2, 'color', [0.3 0.3 0.3], 'MarkerSize', 18, 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',0.1, 'TickLength',[0 0])

%10,000x vs phage
densityplot(result_noidentical_m(:, 3), result_noidentical_m(:,4), [50, 50])
hold on;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
[N, edges, bin] = histcounts(result_noidentical_m(:, 3), nbins);
x = [];
y = [];
sd = [];
for i = 1:1:nbins
    idx = find(bin == i);
    y = [y mean(result_noidentical_m(idx, 4))];
    x = [x mean([edges(i) edges(i+1)])];
    sd = [sd std(result_noidentical_m(idx, 4))];
end
%xlim([0.0013 0.0102])
errorbar(x, y, sd, 'linewidth', 2, 'color', [0.3 0.3 0.3], 'MarkerSize', 18, 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',0.1, 'TickLength',[0 0])

%10,000x vs CAR
densityplot(result_noidentical_m(:, 3), result_noidentical_m(:,6), [50, 50])
hold on;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
[N, edges, bin] = histcounts(result_noidentical_m(:, 3), nbins);
x = [];
y = [];
sd = [];
for i = 1:1:nbins
    idx = find(bin == i);
    y = [y mean(result_noidentical_m(idx, 6))];
    x = [x mean([edges(i) edges(i+1)])];
    sd = [sd std(result_noidentical_m(idx, 6))];
end
%xlim([0.0013 0.0102])
errorbar(x, y, sd, 'linewidth', 2, 'color', [0.3 0.3 0.3], 'MarkerSize', 18, 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',0.1, 'TickLength',[0 0])

%phage vs CAR
densityplot(result_noidentical_m(:, 4), result_noidentical_m(:,6), [50, 50])
hold on;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
[N, edges, bin] = histcounts(result_noidentical_m(:, 4), nbins);
x = [];
y = [];
sd = [];
for i = 1:1:nbins
    idx = find(bin == i);
    y = [y mean(result_noidentical_m(idx, 6))];
    x = [x mean([edges(i) edges(i+1)])];
    sd = [sd std(result_noidentical_m(idx, 6))];
end
%xlim([0.0004 0.01])
errorbar(x, y, sd, 'linewidth', 2, 'color', [0.3 0.3 0.3], 'MarkerSize', 18, 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',0.1, 'TickLength',[0 0])

