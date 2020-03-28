
%3 bio replicate w/ 4 technical replicate each
%missing 60-71 from anderson b/c or corona

outputType = 'gr';

%import files
files = 'no antibiotic data_clean_all.xlsx';

allDat = [];
for i = 1:2
    dat = xlsread(files{i, 1});
    %remove negative values from blanked data
    dat(dat < 0) = 0;
    dat = dat(2:end, :);
    allDat = [allDat; dat];
end 

str = allDat(:, 1);

datFin = allDat(:, 4:end);

s = size(datFin);
%1D median filter (window = 3)
for i = 1:1:s(1)
    datFin(i, :) = medfilt1(datFin(i, :));
end

if(ismember(outputType, {'AUC', 'u+AUC'}))
    % calculate AUC 
    AUC1 = trapz(datFin, 2);
    
    if(ismember(outputType, {'AUC'}))
        datFin = AUC;
    end
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

    
    if(ismember(outputType, {'u', 'u+AUC'}))
        datFin = max(newDat, [], 2);

    end
end

if(ismember(outputType, {'u+AUC'}))
    datFin = [datFin, AUC];
end

%% visualize
biorep = dat(2:13, 2);
bio = unique(biorep);
col = {'r', 'g', 'b'};

str_unique = unique(str);
st = ((length(str_unique)-1) * 3);
en = st+2;
st_e = [[0:3:st]', [2:3:en]'];

figure
for j = 1:length(st_e)
    ind = 0;
    for i = st_e(j, 1):st_e(j, 2)
        ind = ind + 1;
        subplot(2, 8, j)
        hold on;
        plot(1:144, datFin((i*4 + 1):(i*4+4), :), col{ind})
        ylabel(str_unique(j));
    end
end

