%%
% The script is for plotting the growth dynamics of the isolates - the
% clinical library of 203 strains with 4 replicates per strain 
% Figure 2c
%
%%

close all;
clear;

%import dataset with growth rates as the features
data = importClinicalIsolates('gr');

% %indexes of columns corresponding to growth conditions
% one = 2:99;%10,000x
% data_n = normr(data(:, one));
% 
% %margin prediction for 10,000x - 4 fold cv using optimized parameter set
% [saveScores] = singlePredict(data);
% lab = saveScores;
% wrongP = [];
% wrongIdx = [];
% for i = 1:203
%    %if(i ~= lab(i))
%    if(sum(lab(i, :))/4 ~= i)
%        wrongP = [wrongP 1];
%        wrongIdx = [wrongIdx; i];
%    else
%        wrongP = [wrongP 0];
%    end
% end

%Figure 2c - growth rate curves of clinical isolates
figure
hold on;
for i = 0:202
    %subplot(10, 21, i+1)
    start = i*4;
    c = [230 230 248;
        157 150 184;
        119 104 113;
        154 113 151]/250;
    
    subaxis(12, 18,i+1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
    for j = 1:4
        hold on;
        plot(1:98, data_n(start+j, :), 'color', c(j, :), 'linewidth', 0.5);%normalized
        %plot(1:98, data(start+j, :) .*60, 'color', c(j, :), 'linewidth', 1.5);%per hr
    end
    hold on;
%     if(wrongP(i+1) == 1)
%         plot(20, 0.6, 'k*', 'markersize', 5)
%     end
    box off%box on
    axis off 
    axis tight
    ylim([-0.32 0.72])
    set(gca,'linewidth',0.5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])
end

