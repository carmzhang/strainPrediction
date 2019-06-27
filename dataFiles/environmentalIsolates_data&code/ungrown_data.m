%{
This is used to extract all ungrwon strain number
%}

clear;clc;

data = xlsread('ungrown_data.xlsx');
straindata = load('Straindata.mat');

data = data(~isnan(data));
ungrown_strain = unique(sort(data));

%% extract 2 strain number vector of grow and ungrown
AA = ungrown_strain;
BB = straindata.TT_copy;

M = length(BB);
for k = 1:M
    ind = find(AA == BB(k)); % find out data index which has been recorded in straindata
    if(~isempty(ind))
        I(k) = ind;
    else 
        I(k) = 0;
    end
end
%find out the strain index AA contains in BB
I = I(I~=0);
AA(I) = [];%AA: ungrown strain number
save('ungrownstrain.mat','AA');

total = sort([AA;BB]);
number = 1:1:750;

N = length(total);
for m = 1:N
    ind = find(number == total(m));
    if(~isempty(ind))
        II(m) = ind;
    else 
        II(m) = 0;
    end

end
II = II(II~=0);
number(II) = [];
number = number';%number: all number without data

