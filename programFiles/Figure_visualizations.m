%% Visualizations for figures
%
% Miscellaneous plots used to  describe workflow of figures 1, 2, and 4a
%
%%

data = importClinicalIsolates('gr');

%figure 1 growth curve picture
i = [150, 21, 58];
i = 58;
start = i*4 + 1;
c = [230 230 248; 
    157 150 184; 
    119 104 113; 
    154 113 151]/250;
figure
hold on;
for(j = 0:3)
    plot(1:98, data(start+j, 2:99) .*60, 'color', c(j+1, :), 'linewidth', 7);%per hr
end
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

%part a of workflow figure
i = 30;
start = i*4 + 1;
c = [1 1 1; 
    0.8 0.8 0.8; 
    0.6 0.6 0.6; 
    0.2 0.2 0.2];
figure
hold on;
for(j = 0:3)
    plot(1:98, data(start+j, 2:99) .*60, 'color', c(j+1, :), 'linewidth', 7);%per hr
end
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

%part b of workflow figure
i = [30, 103, 194, 57];
c = [157 150 184; 
    124 88 105; 
    154 113 151; 
    11 3 45]/255;

start = i*4 + 1;
for(k = 1:4)
figure
plot(1:98, data(start(k), 2:99) .*60, 'color', c(k, :), 'linewidth', 7);%per hr
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])
end

i_2 = 29;%red
start_2 = i_2*4 + 1;
figure
plot(1:98, data(start_2, 2:99) .*60, 'color', 'r', 'linewidth', 7);%per hr
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

%% part a of correlations explanation
%part a of workflow figure
i = 2;
i_2 = 58;
c = [255 204 51;
    146 209 104;
    119 238 230;
    245 122 202]/255;
start = i*4 + 1;
start_2 = i_2*4 + 1;
figure
plot(1:98, data(start, 2:99) .*60, 'color', c(3, :), 'linewidth', 7);%per hr
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])

figure
plot(1:98, data(start_2, 2:99) .*60, 'color', c(4, :), 'linewidth', 7);%per hr
box on
set(gca,'linewidth',8, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])
