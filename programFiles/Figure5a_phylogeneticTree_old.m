%% Figure 5a plot phylogenetic tree for final figure
% Script plots the phylogenetic tree for 244 clinical isolates and the 
% corresponding antibiotic resistances (true phenotype from 
% clinical microbiology lab testing).
%%

clear;
close all;
load('labels_all.mat')
dat_lbl = [myLabel; myLabelA];

seqs = fastaread('new244bestsnp.fas')
distances = seqpdist(seqs,'method','jukes-cantor','indels','pair');       
phylotree = seqlinkage(distances,'single',seqs)

figure
out = plot(phylotree, 'TerminalLabels', false);
set(gca, 'Visible', 'off')
set(out.BranchLines,'linewidth',2,'Color',[0 0 0])
set(out.BranchDots, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 1)
set(out.axes, 'Color', [1 1 1])
set(out.LeafDots, 'Marker', 'none')
box off
lbel = out.leafNodeLabels;

all_labels = [];
for i = 1:1:244
    
    cur = lbel(i).String;
    cur = erase(cur,'GN');
    cur = erase(cur, 'MDREc');
    cur = str2num(cur);
    all_labels = [all_labels; cur];

end

%% heatmap for resistances
load('resistance_WGSAnalysis_20180501.mat') %with 95% similarity and 50% gene length

comb_CIP = [WGS_CIP4; WGS_CIP4_additional];
SAM = max(WGS_SAM1);
SXT = max(WGS_SXT2);
GM = max(WGS_GM3);
CIP = max(comb_CIP);

threshold = 1;
SAM(SAM < threshold) = 0;
SAM(SAM > threshold) = 1;
SXT(SXT < threshold) = 0;
SXT(SXT > threshold) = 1;
GM(GM < threshold) = 0;
GM(GM > threshold) = 1;
CIP(CIP < threshold) = 0;
CIP(CIP > threshold) = 1;

predictedProfile = [SAM', GM', SXT', CIP'];

all_labels_WGS = [];
for i = 1:1:244
    
    cur = sNames_WGS{i};
    cur = erase(cur,'GN');
    cur = erase(cur, 'MDREc');
    cur = str2num(cur);
    all_labels_WGS = [all_labels_WGS; cur];

end

data = [predictedProfile drugProfile];
data_tree_order = [];
% reorder data to match phylogenetic tree
for i = 1:244
    ind = find(all_labels_WGS == all_labels(i));
    data_tree_order = [data_tree_order; data(ind, :)];
end
n = 2;
R = linspace(238/255, 130/255, n);
G = linspace(221/255, 5/255, n);
B = linspace(255/255, 255/255, n);

%figure 5a
figure
colormap( [R(:), G(:), B(:)] )
imagesc(data_tree_order(:, [5:8]))
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',2,'xticklabel',[], 'yticklabel', [])
