function [roc_all_resfinder, roc_all_card] = resistance_roc_WGS(drug)
%FUNCTION: plots roc curve for wgs analysis for antibiotic of interest
% INPUT:
%   drug - which antibiotic to do prediction for
%          OPTIONS: 1 = SAM   2 = GM   3 = SXT   4 = CIP
%
% OUTPUT: roc curve for wgs for 1 antibiotic
%   roc_all_resfinder - data used to plot roc curve for resfinder
%   roc_all_card - data used to plot roc curve for card

%% roc curve for resfinder
t = readtable('Resfinder_results20190425_fowler and Anderson.xlsx');

vars = {'Sulphonamide' 'Aminoglycoside' 'Beta_lactam' 'Trimethoprim' 'Fluoroquinolone'};
t2 = t{:,vars};
temp = {'none 0 0 |'};
t2(find(strcmp(t2, ''))) = temp;

t{:,vars} = t2;

strains = t{:, 1};

identity_r = zeros(length(t2), 5);
coverage_r = zeros(length(t2), 5);
%split wgs into % identity and coverage
for i = 1:1:5
   for j = 1:length(t2)

       curS = t2{j,i};
       curS = erase(curS, ',');
       cur = strsplit(curS, ' ');
       idx  = find(strcmp(cur, '|'));
       if(length(idx) == 1)
           idx_fin = 2;
       else
           idx_fin = [2 idx(1:end-1)+2];
       end

       %mean of all identities
       cur_identity = cur(idx_fin);
       id = [];
       for k = 1:length(cur_identity)
           id(k) = str2double(cur_identity(k));
       end
       identity_r(j, i) = mean(id);
       
       %mean of all coverages
       cur_coverage = cur(idx-2);
       cov = [];
       for k = 1:length(cur_coverage)
           cov(k) = str2double(cur_coverage(k));
       end
       coverage_r(j, i) = mean(cov);
       
   end
end

load('actual_drug_profile.mat')%SAM, GM, SXT, CIP

%check that strains are same
sNames_p = strrep(sNames_phenotypic, ' ', '');
sNames_w = strrep(strains, ' ', '');
setdiff(sNames_p, sNames_w)
all_strains = [sNames_p sNames_w];
for i = 1:244
    cmp = strcmp(all_strains(i, 1), all_strains(i, 2));
    if(cmp == 0)
       disp(i) 
    end
end

%SAM, GM, SXT, CIP
resfinder_identity = [identity_r(:, [3, 2]), mean(identity_r(:, [1, 4]), 2), identity_r(:, 5)];
resfinder_coverage = [coverage_r(:, [3, 2]), mean(coverage_r(:, [1, 4]), 2), coverage_r(:, 5)];

% threshold_r = [100, 99.9, 99.8, 99.7, 99, 98, 50 0];
% threshold_r2 = [100, 99, 80, 70, 60, 50, 40 30];
threshold_r = [0:1:100];
threshold_r2 = [0:1:100];
roc_all_resfinder = zeros(length(threshold_r), 2);

for i = 1:length(threshold_r)
    temp = resfinder_identity;
    temp(temp < threshold_r(i)) = 0;
    temp(temp >= threshold_r(i)) = 1;
    temp2 = resfinder_coverage;
    temp2(temp2 < threshold_r2(i)) = 0;
    temp2(temp2 >= threshold_r2(i)) = 1;
    
    fin = temp+temp2;
    fin(fin == 1 | fin == 0) = 0;
    fin(fin == 2) = 1; 
    
    predictedProfile = fin;
    
    comb = drugProfile + predictedProfile;
    minus = drugProfile - predictedProfile;
    
    oneTP = sum(comb(:, drug) == 2);
    oneTN = sum(comb(:, drug) == 0);    
    sens = oneTP/sum(drugProfile(:, drug) == 1);
    spec = oneTN/sum(drugProfile(:, drug) == 0);
    roc_all_resfinder(i, :) = [sens, spec];
end

roc_all_resfinder = [roc_all_resfinder;0, 1];
%plot(1-roc_all_resfinder(:, 2), roc_all_resfinder(:, 1), '--o', 'linewidth', 3.5, 'color', 'r', 'markerfacecolor', 'r')

%% roc curve for card
t = readtable('Card_results20190425_fowler and Anderson.xlsx');

vars = {'fluoroquinoloneAntibiotic', 'aminoglycosideAntibiotic', 'sulfonamideAntibiotic', 'TEMBeta_lactamase', 'CMYBeta_lactamase', 'CTX_MBeta_lactamase', 'SHVBeta_lactamase', 'ampC_typeBeta_lactamase'};
t2 = t{:,vars};
t{:,vars} = t2;

strains = t{:, 1};
temp = {'none 0 0 0 |'};
t2(find(strcmp(t2, ''))) = temp;

t{:,vars} = t2;

strains = t{:, 1};
%check that strains are same
sNames_p = strrep(sNames_phenotypic, ' ', '');
sNames_w = strrep(strains, ' ', '');
setdiff(sNames_p, sNames_w)
all_strains = [sNames_p sNames_w];
for i = 1:244
    cmp = strcmp(all_strains(i, 1), all_strains(i, 2));
    if(cmp == 0)
       disp(i) 
    end
end

identity_c = zeros(length(t2), length(vars));
coverage_c = zeros(length(t2), length(vars));
%split wgs into % identity and coverage
for i = 1:1:length(vars)
   for j = 1:length(t2)
       curS = t2{j,i};
       cur = strsplit(curS, ' ');
       idx  = find(strcmp(cur, '|'));
       
       %mean of all identities
       cur_identity = cur(idx-3);
       id = [];
       for k = 1:length(cur_identity)
           id(k) = str2double(cur_identity(k));
       end
       identity_c(j, i) = mean(id);
       
       %mean of all coverages
       cur_coverage = cur(idx-2);
       cov = [];
       for k = 1:length(cur_coverage)
           cov(k) = str2double(cur_coverage(k));
       end
       coverage_c(j, i) = mean(cov);
   end
end

%wgs_prediction_resfinder = [t{:, 5}+t{:, 6}+t{:, 7}+t{:, 8}+t{:, 9} t{:,3}, t{:,4}, t{:,2}];%SAM, GM, SXT, CIP
%SAM, GM, SXT, CIP
card_identity = [mean(identity_c(:, [4, 5, 6, 7, 8]), 2), identity_c(:, [2, 3, 1])];
card_coverage = [mean(coverage_c(:, [4, 5, 6, 7, 8]), 2), coverage_c(:, [2, 3, 1])];

threshold_c = [0:1:100];
threshold_c2 = [0:1:100];
roc_all_card = zeros(length(threshold_c), 2);

for i = 1:length(threshold_c)
    temp = card_identity;
    temp(temp < threshold_c(i)) = 0;
    temp(temp >= threshold_c(i)) = 1;
    temp2 = card_coverage;
    temp2(temp2 < threshold_c2(i)) = 0;
    temp2(temp2 >= threshold_c2(i)) = 1;
    
    fin = temp+temp2;
    fin(fin == 1 | fin == 0) = 0;
    fin(fin == 2) = 1; 
    
    predictedProfile = fin;
    
    comb = drugProfile + predictedProfile;
    minus = drugProfile - predictedProfile;
    
    oneTP = sum(comb(:, drug) == 2);
    oneTN = sum(comb(:, drug) == 0);    
    sens = oneTP/sum(drugProfile(:, drug) == 1);
    spec = oneTN/sum(drugProfile(:, drug) == 0);
    roc_all_card(i, :) = [sens, spec];
end

roc_all_card = [roc_all_card;0, 1];
%plot(1-roc_all_card(:, 2), roc_all_card(:, 1), ':o', 'linewidth', 3.5, 'color', 'r', 'markerfacecolor', 'r')


end

