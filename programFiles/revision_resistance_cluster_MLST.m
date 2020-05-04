%%
% MLST
% - look at sequence type clusters and associated resistance profile
%%

%load tree
distNoRef  = xlsread('Pairwise_distance20191112_MLST.xlsx');

%exp strains for fowler and anderson that are missing
missingF = [2375, 2393, 3253, 3323, 4499, 5443, 5516, 5554, 5992];
missingA = [12, 23, 26, 28, 29];

%remove experimentally missing strains
for i = 1:1:length(missingF)
    distNoRef = distNoRef(distNoRef(:,1)~=missingF(i),:);
    distNoRef = distNoRef(distNoRef(:,2)~=missingF(i),:);
end
for i = 1:1:length(missingA)
    distNoRef = distNoRef(distNoRef(:,1)~=missingA(i),:);
    distNoRef = distNoRef(distNoRef(:,2)~=missingA(i),:);
end

%load strain name with resistance profile
%1. label 
%2. SAM
%3. GM
%4. SXT
%5. CIP
resistance_label  = xlsread('Strains_label_resistance.xlsx');
isolates = unique([resistance_label(:, 1)]);

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

%for clustered strains compile resistances
cluster_resistance = cell(s(1), 1);
for i = 1:1:s(1)
    cur = clusteredStrains{i, 1};
    cluster_resistance{i, 1} = [];
    for j = 1:length(cur)
       temp = cur(j); 
       idx = find(resistance_label(:, 1) == temp);
       lab = resistance_label(idx,:);
       cluster_resistance{i, 1} = [cluster_resistance{i, 1}; lab];
    end
end

allClusters = [];
for i = 1:1:s(1)
    allClusters = [allClusters; cluster_resistance{i, 1}; 2 2 2 2 2];
end


map = [1 1 1;
       0 0 1;
       0 0 0];
HeatMap(allClusters(:, 2:end), 'Colormap', 'redbluecmap')


%ab response for PIPTAZ, TIM, SAM
%1. PIPTAZ - piperacillin-tazobactam
%2. TIM - ticarcillin-clavulanate
%3. SAM - ampicillin-sulbactam

load('fowler_ab_supp_revision.mat')
%HeatMap(AbresponsefowlerisolatesS2, 'Colormap', 'redbluecmap')
str = '#E0765F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str2 = '#a20021';
color2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
map = [0 0 0;
       1 1 1;
       color];
imagesc(AbresponsefowlerisolatesS2)
colormap(map)
axis off;

%compare PIPTAZ and SAM
temp = AbresponsefowlerisolatesS2(11:end, 1) == AbresponsefowlerisolatesS2(11:end, 3);
mean(temp)

%compare TIM and SAM
temp2 = AbresponsefowlerisolatesS2(1:183, 2) == AbresponsefowlerisolatesS2(1:183, 3);
mean(temp2)

%compare PIPTAZ and TIM
temp3 = AbresponsefowlerisolatesS2(11:183, 1) == AbresponsefowlerisolatesS2(11:183, 2);
mean(temp3)

