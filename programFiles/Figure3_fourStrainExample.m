%% Plot the 4 strain example for figure 3
% The script generates the plots for Figure 3a-b, which is a 4 strain 
% example of the svm training process. Here, we extract the top 2 
% components from pca to use as the features.
% (1-4) growth rate curve for the 4 strains (Figure 3a)
% (5) plot full dataset for 4 strain example
% (6-9) plot svm boundary for each of the 4 classes (Figure 3b)
% (10) plot training set 
%%

data = importClinicalIsolates('gr');

i = [15, 22, 202, 200];
start = i*4 + 1;
fin = i*4+4;

%color code strains
c = [255 204 51;
    146 209 104;
    119 238 230;
    245 122 202]/255;
for(j = 1:4)
    figure;
    plot(1:98, data(start(j):fin(j), 2:99) .*60, 'color', c(j, :), 'linewidth', 2);%per hr
    box on
    ylim([0, 0.3])
    if(j > 1)
        set(gca,'linewidth',5, 'xticklabel', [], 'xtick', [], 'yticklabel', [], 'ytick', [])
    else
        set(gca,'fontsize', 36, 'FontWeight','bold', 'linewidth',5)
    end
end

%%%%%%%%%

%plot PCA of 4 examples
one = 2:99;
i = [15, 22, 202, 200];
start = i*4 + 1;
fin = i*4+4;
[coeff score latent] = pca(data([start(1):fin(1), start(2):fin(2), start(3):fin(3), start(4):fin(4)], one).*60);
figure
hold on
for i = 1:16
    if(i <=4)
        scatter(score(i, 1), score(i, 2), 100, [255 204 51]/255, 'filled')
    elseif(i>4 && i<=8)
        scatter(score(i, 1), score(i, 2), 100, [146 209 104]/255, 'filled')
    elseif(i>8 && i<=12)
        scatter(score(i, 1), score(i, 2), 100, [119 238 230]/255, 'filled')
    else
        scatter(score(i, 1), score(i, 2), 100, [245 122 202]/255, 'filled')
    end
end
set(gca, 'fontsize', 26, 'FontWeight','bold','linewidth',4)
hold off;

%subset - train
s = size(score);
i = 1;

temp = 1:1:s(1);
temp(i:4:end) = [];

trainX = score(temp, 1:2);
testX = score(i:4:end, 1:2);

N_test = 4;% number of testing samples
N_train = 4*3;% number of training samples
K = 4;%number of classes
M = 2;%dimension of features

test_label = [1 2 3 4]';
train_label = [1 1 1 2 2 2 3 3 3 4 4 4]';
test_feature = testX;
train_feature = trainX;

%run svm model
[predict_label,result]=multi_class_svm_plot(test_label,train_label,test_feature,train_feature,N_test,N_train, 'rbf', 1, 0, 1);

%plot training set
figure;
hold on;
for j = 1:12
    if(j <=3)
        col = [255 204 51]/255;
    elseif(j>3 && j<=6)
        col = [146 209 104]/255;
    elseif(j>6 && j<=9)
        col = [119 238 230]/255;
    else
        col = [245 122 202]/255;
    end
    scatter(trainX(j, 1), trainX(j, 2), 600, 'MarkerFaceColor', col, 'MarkerEdgeColor', [1 1 1])
end
ylim([-0.1 0.15])
xlim([-0.25 0.25])
set(gca, 'fontsize', 26, 'FontWeight','bold','linewidth',4)
legend('off')
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
hold off;
