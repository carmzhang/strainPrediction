%%
% Plot phylogenetic trees of environmental isolates - by order
% 
%%

%% CHOOSE ORDER OF INTEREST HERE
% OPTIONS: 1 - Actinomycetales
%          2 - Bacillales [most isolates in this order]
%          3 - Burkholderiales
%          4 - Corynebacteriales
%          5 - Enterobacteriales
%          6 - Eukaryota
%          7 - Lactobacillales
%          8 - Micrococcales
%          9 - Neisseriales
%          10 - Pseudomonadales
%          11 - Rhizobiales
%          12 - Streptomycetales
%          13 - unclassified
%          14 - Xanthomonadales
%          15 - reference
%%
%Ex: tree below corresponds to Actinomycetales
filenames = {'Actinomycetales_bestsnp.fasta';
              'Bacillales_bestsnp.fasta';
              'Burkholderiales_bestsnp.fasta';
              'Corynebacteriales_bestsnp.fasta';
              'Enterobacteriales_bestsnp.fasta';
              'Eukaryota_bestsnp.fasta';
              'Lactobacillales_bestsnp.fasta';
              'Micrococcales_bestsnp.fasta';
              'Neisseriales_bestsnp.fasta';
              'Pseudomonadales_bestsnp.fasta';
              'Rhizobiales_bestsnp.fasta';
              'Streptomycetales_bestsnp.fasta';
              'unclassified_bestsnp.fasta';
              'Xanthomonadales_bestsnp.fasta';
              'reference_20190310.fasta'};
file = filenames{15};%change index here to change order

%import phylogenetic tree
seqs = fastaread(file)
distances = seqpdist(seqs,'method','jukes-cantor','indels','pair');       
phylotree = seqlinkage(distances,'single',seqs)

figure
%out = plot(phylotree, 'TerminalLabels', false);
out = plot(phylotree);
set(gca, 'Visible', 'off')
set(out.BranchLines,'linewidth',2,'Color',[0 0 0])
set(out.BranchDots, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 1)
set(out.axes, 'Color', [1 1 1])
set(out.LeafDots, 'Marker', 'none')
box off
lbel = out.leafNodeLabels;

%all_labels are the FS# for the isolates in the order that they appear from
%top to bottom  in the phylogenetic tree.
all_labels = [];
for i = 1:1:length(lbel)
    
    cur = lbel(i).String;
    cur = erase(cur,'FS');
    cur = str2num(cur);
    all_labels = [all_labels; cur];

end
