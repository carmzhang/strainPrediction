%%
% The script is for plotting the growth dynamics of the isolates - the
% environmental isolate library of 143 strains with 12 replicates per
% strain.
% Supplementary Figure 2.2
%
%%

close all;
clear;

%% import environmental isolates
data = importEnvironmentalIsolates(1);%import data with new labels
data_n = normr(data(:, 2:end));

%Supplementary Figure 2.2 - growth rate curves of environmental isolates
figure
hold on;
for i = 0:142
    start = i*12;
    c = [linspace(230, 154, 12)', linspace(230, 113, 12)', linspace(248, 151, 12)';
       ]/255;
    
    subaxis(15, 10,i+1, 'Spacing', 0.02, 'Padding', 0, 'Margin', 0);
    for j = 1:12
        hold on;
        plot(1:98, data_n(start+j, :), 'color', c(j, :), 'linewidth', .5);%normalized
    end
    box on
    set(gca,'linewidth',0.5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

    ylim([-0.7 0.91])
    box off%box on
    axis off 
    axis tight
    set(gca,'linewidth',0.5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])


end

