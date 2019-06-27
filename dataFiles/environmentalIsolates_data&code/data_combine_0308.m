%{
used to combine all data: 12, 18, 24 replicates
%}

clear;clc;

data1 = load('Straindata_new_0308new.mat');
data2 = load('Straindata_trip_0308.mat');
data3 = load('Straindata_multi_0308.mat');

strain_rep = data1.Data_copy;
strain_trip = data2.Straindata_trip;
strain_multi = data3.Straindata_multi;

trip_num = size(strain_trip,1)/18;
multi_num = size(strain_multi,1)/24;

for i = 1:trip_num
    AA = strain_trip(1+(i-1)*18:i*18-6,:);
    strain_rep = [strain_rep;AA];
end

for k = 1:multi_num
    BB = strain_multi(1+(k-1)*24:k*24-12,:);
    strain_rep = [strain_rep;BB];
end

strain_rep = sortrows(strain_rep);
Data_copy = strain_rep;


%delete the extra 767 768
ind_1 = find(Data_copy(:,1) == 767);
ind_2 = find(Data_copy(:,1) == 768);
Data_copy([ind_1,ind_2],:) = [];
TT_copy = unique(Data_copy(:,1));

%% save the latest file
save('Straindata_latest_0308.mat','Data_copy','TT_copy');
%% plot 12 replicates
for k = 1:100
    figure(1);
    subplot(10,10,k);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end
% 
for k = 101:200
    figure(2);
    subplot(10,10,k-100);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end 

for k = 201:300
    figure(3);
    subplot(10,10,k-200);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end 
for k = 301:400
    figure(4);
    subplot(10,10,k-300);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end 
for k = 401:500
    figure(5);
    subplot(10,10,k-400);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end 
for k = 501:600
    figure(6);
    subplot(10,10,k-500);
    plot(Data_copy(1+(k-1)*12:k*12-6,2:100)','k');
    hold on;
    plot(Data_copy(k*12-5:k*12,2:100)','r');
    title(TT_copy(k));
    xticks([]);
    ylim([0 0.8]);
end 



