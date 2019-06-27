function [data] = importEnvironmentalIsolates(lab)
%FUNCTION: returns the dataset for environmental isolates with all 
%          unique strains 
%          *features are growth rates
% INPUT:
%   lab - what the first column (label) will look like
%           0 for original bdr# labels
%           1 for 1 -> last label (used for prediction)
% OUTPUT: data - dataset of all unique strains with 4 replicates per strain


ord = readtable('Order_list_seq_20190312.xlsx', 'TreatAsEmpty', 'N/A');
fw = 199:229;%different set of isolates - remove

o_seq = cell(width(ord), 1);
for i = 1:1:width(ord)
    o1 = ord{:, i};
    o1(isnan(o1)) = [];
    o_seq{i} = o1;
end

o_all = [];
for i = 1:1:width(ord)
    o_all = [o_all; o_seq{i}];
end

%check for duplicates
A = o_all;
[uniqueA i j] = unique(A,'first');
indexToDupes = find(not(ismember(1:numel(A),i)))

unique_isolates_order = [15, 18, 7, 586, 154, 335, 60, 157, 74, 363, 367];%FS#

%GC data
load('Straindata_latest_0308.mat') %650 isolates
%TT_copy are the unique isolate # (bdr#)
%Data_copy is GC in replicates of 12

%match ID of GC with FSID (sequence ID)
match_id = readtable('Environmental isolates_information_NEW.xlsx', 'TreatAsEmpty', 'N/A');
match_fin = [match_id{:, 2}, match_id{:, 5}];%column 1: bdr#, 2: FS#

%convert FS to bdr#: o_all -> o_bdr
o_all = sort(o_all);
o_bdr = [];
for i = 1:1:length(o_all)
    idx = find(match_fin(:, 2) == o_all(i));
    o_bdr = [o_bdr; match_fin(idx, 1)];
end
total_isolates = intersect(o_bdr, TT_copy);

%% ------------cluster identical isolates

[~,sheet_name]=xlsfinfo('all_results_20190312_mod.xlsx');
sheet_name = sheet_name(2:end-2);
data = cell(numel(sheet_name), 1);
for k = 1:numel(sheet_name)
    data{k} = xlsread('all_results_20190312_mod.xlsx',sheet_name{k});
end

%cluster identical strains - distance = 0
for j = 1:numel(sheet_name)
    cur = data{j};
    identical = cur(find(cur(:, 3) == 0), :);
    
    s_i = size(identical);
    identical_bdr = [];
    for i = 1:1:s_i(1)
        idx_1 = find(match_fin(:, 2) == identical(i, 1));
        idx_2 = find(match_fin(:, 2) == identical(i, 2));
        identical_bdr = [identical_bdr; match_fin(idx_1, 1), match_fin(idx_2, 1)];
    end
    identical_bdr_all{j} = identical_bdr;
    
    if(~isempty(identical_bdr))
        strains = cell(s_i(1), 1);
        found = [];
        prev = identical_bdr(end, 2);
        count = 1;
        for i = 0:1:s_i(1)-1
            temp = identical_bdr(end-i, 1:2);
            if(temp(1, 2) == prev || ~isempty(find(strains{count, 1} == temp(1, 2))) )
                strains{count, 1} = unique([strains{count, 1}, temp]);
                found = [found, temp];
            else
                count = count + 1;
                strains{count, 1} = unique([strains{count, 1}, temp]);
            end
            prev = temp(1, 2);
        end
        strains = strains(~cellfun('isempty',strains));
        
        %strFin - each cell represents a strain cluster of isolates
        strFin = strains;
        for i_in = 1:1:numel(strFin)
            for j_in = i_in+1:1:numel(strFin)
                if(j_in ~= i_in)
                    if( ~isempty(intersect(strains{i_in}, strains{j_in})) && ~isempty(strFin{i_in}))%share elements
                        strFin{i_in, 1} = unique([strains{i_in}, strains{j_in}]);
                        strFin{j_in, 1} = [];
                        
                    end
                    
                end
            end
        end
        strFin = strFin(~cellfun('isempty',strFin));
        
        strains_all{j} = strFin;
    else
        strains_all{j} = [];
    end
    
end
identical_bdr_dat = identical_bdr_all(~cellfun('isempty',identical_bdr_all));
strains_dat = strains_all(~cellfun('isempty',strains_all));

%extract unique strains from clusters
%o_bdr - all isolates
included_str = [];%counting isolates in clusters
with_gc = [];%isolates in clusters with corresponding gc
need_dat  = [];%isolates in clusters w/o gc
for i = 1:numel(strains_dat)
    
    cur =  strains_dat{i};
    for j = 1:numel(cur)
        cur_in  = cur{j};
        included_str = [included_str, cur_in];
        idx  = intersect(cur_in, TT_copy);
        
        if(isempty(idx))
            need_dat = [need_dat; cur_in];
        else
            with_gc = [with_gc; idx(1)];
        end
    end
    
end

nonclustered = setdiff(o_bdr, included_str);
nonclustered_gc = intersect(nonclustered, TT_copy);

%% need data:
%setdiff(nonclustered, nonclustered_gc)
%%%%%%%%%%%%%

data_str = sort([nonclustered_gc; with_gc]);

[tf loc] = ismember(Data_copy(:, 1), data_str);
idx = find(tf);
data_seq = Data_copy(idx, 2:end);
data_seq_lab = Data_copy(idx, 1);
%-----------------------------------------------------------
%convert raw data to growth rates
s = size(data_seq);
datFin = [];
%1D moving average filter (window = 5)
for i = 1:1:s(1)
    data_seq(i, :) = conv2(data_seq(i, :), ones(1, 5)/5, 'same');
end
data_seq = data_seq - min(min(data_seq))+0.0000001;

%convert raw data to derivative
time = 10/60;
data_seq_tr = gradient(log(data_seq), time);
s = size(data_seq_tr);

%randomize of order of replicates
data_FIN = [];
for i = 1:1:s(1)/12
    s_in = RandStream('mt19937ar','Seed',i);
    idx = randperm(s_in, 12);%12 replicates per isolate
    start = (i-1).*12 + 1;
    idx_fin = (idx-1) + start;
    temp = data_seq_tr(idx_fin,:);
    data_FIN = [data_FIN; temp];
end

v = 1:length(data_FIN)/12;
if(lab == 0)

    data = [data_seq_lab, data_FIN];
else
    labels = repelem(v,12)'; %for prediction ONLY
    data = [labels, data_FIN];
end


end

