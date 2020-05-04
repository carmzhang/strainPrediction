%%
% The script is for plotting the growth dynamics of the isolates - the
% environmental isolate library of 143 strains with 12 replicates per
% strain.
% Supplementary Figure 2.2
%
%%

close all;
clear;

%% import clinical isolates
outputType = 'gr';
seq = 'snp';%other option is mlst

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
%1D median filter (window = 3)
for i = 1:1:s(1)
    datFin(i, :) = medfilt1(datFin(i, :));
end

if(ismember(outputType, {'gr', 'u', 'u+AUC'}))
    %alternative: convert raw data to derivative
    time = 0:10:(s(2)-1)*10;
    newDat1 = [];

    for i = 1:1:s(1)
        x1 = time(1:end-1)';
        x2 = time(2:end)';
        y1 = datFin(i, 1:end-1)';
        y2 = datFin(i, 2:end)';
        slopes = (y2 - y1) ./ (x2 - x1);
        newDat(i, :) = slopes';
       
    end
    datFin = newDat;
end
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
%data_n = normr(data(:, 2:end));
data_n = data;
%Supplementary Figure 1.6 - growth rate curves of clinical isolates - snp
%strains
figure
hold on;
for i = 0:193
    start = i*12;

    subaxis(15, 13,i+1, 'Spacing', 0.02, 'Padding', 0, 'Margin', 0);
    for j = 1:12
        hold on;
        if(j>= 1 && j <= 4)
            c = 'r';
        elseif(j>= 5 && j <= 8)
            c = 'g';
        else
            c = 'b';
        end
        plot(1:144, data_n(start+j, 2:end), 'color', c, 'linewidth', .8);%normalized
    end
    box on
    set(gca,'linewidth',0.5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

    %ylim([-0.7 0.91])
    box off%box on
    axis off 
    axis tight
    set(gca,'linewidth',0.5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

end

