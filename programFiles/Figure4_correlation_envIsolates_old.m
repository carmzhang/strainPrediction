%% 
% The  script plots the right subplot in Figure 4b, the correlation between 
% phenotype and genotype for the environmental isolates. Here, we also 
% display the spearman coefficient and corresponding p-value (significance)
% and heat map for all orders found in the library of environmental
% isolates. The order you want to see these values for can be chosen in line
% 121. Figure 4, Supplementary Figure 2.4
%
%% environmental isolates 
clear;
close all;

%import sequence information - paired distances
[~,sheet_name]=xlsfinfo('all_results_20190312_mod.xlsx');
sheet_name = sheet_name(2:end-2);
data_pairDist = cell(numel(sheet_name), 1);
isolates_perOrder = cell(numel(sheet_name), 1);
for k = 1:numel(sheet_name)
    temp = xlsread('all_results_20190312_mod.xlsx',sheet_name{k});
    isolates_perOrder{k} = unique([temp(:, 1), temp(:, 2)]);
    data_pairDist{k} = temp;
end

%GC data
load('Straindata_latest_0308.mat') %GC for 650 isolates
%TT_copy are the unique isolate # (bdr#)
%Data_copy is GC in replicates of 12
data = Data_copy;
data_m = reshape(mean(reshape(data,12,[])),[],100);%mean of 12 replicates

byOrder_wgs_gc = cell(numel(sheet_name), 2);
for k = 1:numel(sheet_name)
    
    %per order - correlation between genotype and phenotype
    isolates = isolates_perOrder{k};
    isolates = intersect(isolates, TT_copy);
    seq_pair = data_pairDist{k};
    
    pair_wgs = zeros(length(isolates), length(isolates));
    pair_phenotype_gc = cell(length(isolates), length(isolates));
    pair_phenotype_m = zeros(length(isolates), length(isolates));
    
    for i = 1:length(isolates)
        
        i1 = isolates(i);
        
        for j = 1:length(isolates)
            i2 = isolates(j);
            if(i1 ~= i2)
                idx1_m = find(data_m(:, 1) == i1);
                idx2_m = find(data_m(:, 1) == i2);
                pair_phenotype_m(i, j) = pdist([data_m(idx1_m, 2:end); data_m(idx2_m, 2:end)]);
                pair_phenotype_gc{i, j} = [data_m(idx1_m, 2:end); data_m(idx2_m, 2:end)];
                
                if(~isempty(intersect(find(seq_pair(:, 1) == i1),find(seq_pair(:, 2) == i2))))
                    row_in = seq_pair(intersect(find(seq_pair(:, 1) == i1),find(seq_pair(:, 2) == i2)), 3);
                    pair_wgs(i, j) = row_in;
                else
                    row_in = seq_pair(intersect(find(seq_pair(:, 2) == i1),find(seq_pair(:, 1) == i2)), 3);
                    pair_wgs(i, j) = row_in;
                end
            end
            
        end
    end    
    byOrder_wgs_gc{k, 2} = pair_phenotype_gc; 
    byOrder_wgs_gc{k, 3} = pair_wgs;
    byOrder_wgs_gc{k, 4} = pair_phenotype_m;

    result_noidentical_m = [];
    for i = 1:length(isolates)
        for j = 1:length(isolates)
            if(i ~= j)
                if(pair_wgs(i, j) ~= -214748364)
                    result_noidentical_m = [result_noidentical_m; pair_wgs(i, j), pair_phenotype_m(i, j)];
                end
            end
        end
    end
    byOrder_wgs_gc{k, 1} = result_noidentical_m;
end

%% Supplementary Figure 2.4
% Density plots for all orders of environmental isolate library
for i = 1:numel(sheet_name)
    
    if(~isempty(byOrder_wgs_gc{i, 1}))
        densityplot(byOrder_wgs_gc{i, 1}(:, 1), byOrder_wgs_gc{i, 1}(:, 2), [50, 50])%wgs vs 10,000x
        ax = gca;
        ax.XRuler.Exponent = 0;
        ax.YRuler.Exponent = 0;
        title(sheet_name{i})
        set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)
    end
end

%% Figure 4 
% -  mean of 12 replicates density plot
nbins = 8;
%for figure 3 (bacillales)
densityplot(byOrder_wgs_gc{4, 1}(:, 1), byOrder_wgs_gc{4, 1}(:, 2), [50, 50])%wgs vs 10,000x
hold on;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0; 
title(sheet_name{4})
[N, edges, bin] = histcounts(byOrder_wgs_gc{4, 1}(:, 1), nbins);
x = [];
y = [];
sd = [];
for i = 1:1:nbins
    idx = find(bin == i);
    y = [y mean(byOrder_wgs_gc{4, 1}(idx, 2))];
    x = [x mean([edges(i) edges(i+1)])];
    sd = [sd std(byOrder_wgs_gc{4, 1}(idx, 2))];
end
errorbar(x, y, sd, 'linewidth', 2, 'color', [0.3 0.3 0.3], 'MarkerSize', 18, 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0])
xlim([0.015 1.55])
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',0.1, 'TickLength',[0 0])

%% CHOOSE order to look at correlation between phylogeny and phenotype 
% OPTIONS: 1 - Actinomycetales
%          3 - Micrococcales
%          4 - Bacillales
%          5 - Lactobacillales
%          7 - Streptomycetales
%          9 - Enterobacteriales
%          10 - Pseudomonadales

%Ex: chosen is oder 10, Pseudomonadales
order = 10;%change this line to change order
[out pval] = bramila_mantel(byOrder_wgs_gc{order, 3},byOrder_wgs_gc{order, 4},5000,'spearman')



