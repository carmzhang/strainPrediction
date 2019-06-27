%% Figure 5a
% The script plots with lines through circle image of Figure 5a.
%%

clear;
close all;

%% for resistances - condition combination 1->4
%1:accuracy 2: trp 3: tnr
load('resistance_predictions_FINAL_20181018.mat')
CIP = CIP([1:4, 6:end, 5], :);
GM = GM([1:4, 6:end, 5], :);
SAM = SAM([1:4, 6:end, 5], :);
SXT = SXT([1:4, 6:end, 5], :);

%using collection of genes
SAM_wgs = [70.49, 93.33, 42.20];%acc, tpr, tnr
GM_wgs = [55.74, 100, 47.32];
SXT_wgs = [86.07, 92.14, 77.88];
CIP_wgs = [88.52, 83.56, 95.92];

%using Resfinder
SAM_wgs_1 = [69.67, 88.89, 45.87];
GM_wgs_1 = [55.33, 94.87 47.80];
SXT_wgs_1 = [86.48, 92.14, 78.85];
CIP_wgs_1 = [53.69 22.6 100];

%using CARD
SAM_wgs_2 = [69.26, 89.63, 44.04];
GM_wgs_2 = [17.21, 100, 1.46];
SXT_wgs_2 = [85.25, 83.57, 87.50];
CIP_wgs_2 = [59.84, 100, 0];

x0 = 0;
y0 = 0;
%top right quadrant - SAM
angle_SAM = linspace(90,0, 49);
angle_SAM(find(angle_SAM == 0)) = [];
angle_SAM(find(angle_SAM == 30)) = [];
angle_SAM(find(angle_SAM == 60)) = [];
angle_SAM(find(angle_SAM == 90)) = [];

%label growth conditions - tnr for SAM
cond_label = fliplr(angle_SAM(end-14:end));
x_cond_1 = cosd(cond_label([1, 3, 4, 5, 9, 10, 11, 15]))' *1.05;%10,000x
y_cond_1 = sind(cond_label([1, 3, 4, 5, 9, 10, 11, 15]))' *1.05;

x_cond_2 = cosd(cond_label([1, 2, 4, 5, 7, 8, 11, 14]))' *1.12;%phage
y_cond_2 = sind(cond_label([1, 2, 4, 5, 7, 8, 11, 14]))' *1.12;

x_cond_3 = cosd(cond_label([1, 2, 3, 5, 6 8, 10, 13]))' *1.19;%100x
y_cond_3 = sind(cond_label([1, 2, 3, 5, 6 8, 10, 13]))' *1.19;

x_cond_4 = cosd(cond_label([1, 2, 3, 4, 6 7, 9, 12]))' *1.26;%carb
y_cond_4 = sind(cond_label([1, 2, 3, 4, 6 7, 9, 12]))' *1.26;

SAM_reshape = [SAM(:, 1); SAM(:, 2); SAM(:, 3)]/100;
x_SAM = cosd(angle_SAM)'.*SAM_reshape;
y_SAM = sind(angle_SAM)'.*SAM_reshape;

x_SAM_tnr = cosd(0:0.1:30)'.*[SAM_wgs(3), SAM_wgs_1(3), SAM_wgs_2(3)]/100;
y_SAM_tnr = sind(0:0.1:30)'.*[SAM_wgs(3), SAM_wgs_1(3), SAM_wgs_2(3)]/100;
x_SAM_tpr = cosd(30:0.1:60)'.*[SAM_wgs(2), SAM_wgs_1(2), SAM_wgs_2(2)]/100;
y_SAM_tpr = sind(30:0.1:60)'.*[SAM_wgs(2), SAM_wgs_1(2), SAM_wgs_2(2)]/100;
x_SAM_acc = cosd(60:0.1:90)'.*[SAM_wgs(1), SAM_wgs_1(1), SAM_wgs_2(1)]/100;
y_SAM_acc = sind(60:0.1:90)'.*[SAM_wgs(1), SAM_wgs_1(1), SAM_wgs_2(1)]/100;

%top left quadrant - CIP
angle_CIP = linspace(180,90, 49);
angle_CIP(find(angle_CIP == 90)) = [];
angle_CIP(find(angle_CIP == 120)) = [];
angle_CIP(find(angle_CIP == 150)) = [];
angle_CIP(find(angle_CIP == 180)) = [];

CIP_reshape = [CIP(:, 1); CIP(:, 2); CIP(:, 3)]/100;
x_CIP = cosd(angle_CIP)'.*CIP_reshape;
y_CIP = sind(angle_CIP)'.*CIP_reshape;

x_CIP_tnr = cosd(90:0.1:120)'.*[CIP_wgs(3), CIP_wgs_1(3), CIP_wgs_2(3)]/100;
y_CIP_tnr = sind(90:0.1:120)'.*[CIP_wgs(3), CIP_wgs_1(3), CIP_wgs_2(3)]/100;
x_CIP_tpr = cosd(120:0.1:150)'.*[CIP_wgs(2), CIP_wgs_1(2), CIP_wgs_2(2)]/100;
y_CIP_tpr = sind(120:0.1:150)'.*[CIP_wgs(2), CIP_wgs_1(2), CIP_wgs_2(2)]/100;
x_CIP_acc = cosd(150:0.1:180)'.*[CIP_wgs(1), CIP_wgs_1(1), CIP_wgs_2(1)]/100;
y_CIP_acc = sind(150:0.1:180)'.*[CIP_wgs(1), CIP_wgs_1(1), CIP_wgs_2(1)]/100;

%bottom left quadrant - SXT
angle_SXT = linspace(270,180, 49);
angle_SXT(find(angle_SXT == 180)) = [];
angle_SXT(find(angle_SXT == 210)) = [];
angle_SXT(find(angle_SXT == 240)) = [];
angle_SXT(find(angle_SXT == 270)) = [];

SXT_reshape = [SXT(:, 1); SXT(:, 2); SXT(:, 3)]/100;
x_SXT = cosd(angle_SXT)'.*SXT_reshape;
y_SXT = sind(angle_SXT)'.*SXT_reshape;

x_SXT_tnr = cosd(180:0.1:210)'.*[SXT_wgs(3), SXT_wgs_1(3), SXT_wgs_2(3)]/100;
y_SXT_tnr = sind(180:0.1:210)'.*[SXT_wgs(3), SXT_wgs_1(3), SXT_wgs_2(3)]/100;
x_SXT_tpr = cosd(210:0.1:240)'.*[SXT_wgs(2), SXT_wgs_1(2), SXT_wgs_2(2)]/100;
y_SXT_tpr = sind(210:0.1:240)'.*[SXT_wgs(2), SXT_wgs_1(2), SXT_wgs_2(2)]/100;
x_SXT_acc = cosd(240:0.1:270)'.*[SXT_wgs(1), SXT_wgs_1(1), SXT_wgs_2(1)]/100;
y_SXT_acc = sind(240:0.1:270)'.*[SXT_wgs(1), SXT_wgs_1(1), SXT_wgs_2(1)]/100;

%bottom right quadrant - GM
angle_GM = linspace(360,270, 49);
angle_GM(find(angle_GM == 270)) = [];
angle_GM(find(angle_GM == 300)) = [];
angle_GM(find(angle_GM == 330)) = [];
angle_GM(find(angle_GM == 360)) = [];

GM_reshape = [GM(:, 1); GM(:, 2); GM(:, 3)]/100;
x_GM = cosd(angle_GM)'.*GM_reshape;
y_GM = sind(angle_GM)'.*GM_reshape;

x_GM_tnr = cosd(270:0.1:300)'.*[GM_wgs(3), GM_wgs_1(3), GM_wgs_2(3)]/100;
y_GM_tnr = sind(270:0.1:300)'.*[GM_wgs(3), GM_wgs_1(3), GM_wgs_2(3)]/100;
x_GM_tpr = cosd(300:0.1:330)'.*[GM_wgs(2), GM_wgs_1(2), GM_wgs_2(2)]/100;
y_GM_tpr = sind(300:0.1:330)'.*[GM_wgs(2), GM_wgs_1(2), GM_wgs_2(2)]/100;
x_GM_acc = cosd(330:0.1:360)'.*[GM_wgs(1), GM_wgs_1(1), GM_wgs_2(1)]/100;
y_GM_acc = sind(330:0.1:360)'.*[GM_wgs(1), GM_wgs_1(1), GM_wgs_2(1)]/100;

%best prediction
[val idx_SAM] = max(mean([SAM(:, 1), SAM(:, 2), SAM(:, 3)], 2));
[val idx_GM] = max(mean([GM(:, 2)], 2));
[val idx_SXT] = max(mean([SXT(:, 1), SXT(:, 2), SXT(:, 3)], 2));
[val idx_CIP] = max(mean([CIP(:, 1), CIP(:, 2), CIP(:, 3)], 2));

%color 
a = [22 135 0]/255;%tnr
b = [0 97 255]/255;%tpr
c = [178 53 0]/255;%acc
d = [0 0 0;
    0 0 0;
    0 0 0];%all same colors
line_t = {'-', '--', ':'};

wid_arc = 1.5;

%% start figure
figure
hold on;
%test condition lines
tnr_a = angle_SAM(end-14:end);
x = cosd(tnr_a)'.*1.3;
y = sind(tnr_a)'.*1.3;
for i = 1:length(tnr_a)
    h = plot([0 x(i)], [0 y(i)], 'linewidth', 0.2, 'Color', [0 0 0], 'LineStyle', ':');
end

%condition label
scatter(x_cond_1, y_cond_1, 50, [114, 255, 145]/255, 'filled')
scatter(x_cond_2, y_cond_2, 50, [149, 226, 245]/255, 'filled')
scatter(x_cond_3, y_cond_3, 50, [216, 195, 251]/255, 'filled')
scatter(x_cond_4, y_cond_4, 50, [243, 202, 203]/255, 'filled')

for i = 1:length(x_SAM)
    if(i <= 15)
        co = [252 159 91]/255;
        if(i == idx_SAM)
            w = 4;
            scatter(cosd(angle_SAM(i))*(SAM_reshape(i)+0.03), sind(angle_SAM(i))*(SAM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    elseif(i >15 && i <= 30)
        co = [86 238 244]/255;
        if(i == idx_SAM+15)
            w = 4;
            scatter(cosd(angle_SAM(i))*(SAM_reshape(i)+0.03), sind(angle_SAM(i))*(SAM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    else
        co = [50 232 117]/255;
        if(i == idx_SAM+30)
            w = 4;
            scatter(cosd(angle_SAM(i))*(SAM_reshape(i)+0.03), sind(angle_SAM(i))*(SAM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    end
     h = plot([0 x_SAM(i)], [0 y_SAM(i)], 'linewidth', w, 'Color', co);
end

for k1 = 1:3
    plot(x_SAM_tnr(:, k1), y_SAM_tnr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_SAM_tpr(:, k1), y_SAM_tpr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_SAM_acc(:, k1), y_SAM_acc(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
end

for i = 1:length(x_CIP)
    if(i <= 15)
        co = [252 159 91]/255;
        if(i == idx_CIP)
            w = 4;
            scatter(cosd(angle_CIP(i))*(CIP_reshape(i)+0.03), sind(angle_CIP(i))*(CIP_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    elseif(i >15 && i <= 30)
        co = [86 238 244]/255;
        if(i == idx_CIP+15)
            w = 4;
            scatter(cosd(angle_CIP(i))*(CIP_reshape(i)+0.03), sind(angle_CIP(i))*(CIP_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    else
        co = [50 232 117]/255;
        if(i == idx_CIP+30)
            w = 4;
            scatter(cosd(angle_CIP(i))*(CIP_reshape(i)+0.03), sind(angle_CIP(i))*(CIP_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    end
     h = plot([0 x_CIP(i)], [0 y_CIP(i)], 'linewidth', w, 'Color', co);
end

for k1 = 1:3
    plot(x_CIP_tnr(:, k1), y_CIP_tnr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_CIP_tpr(:, k1), y_CIP_tpr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_CIP_acc(:, k1), y_CIP_acc(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
end

for i = 1:length(x_SXT)
    if(i <= 15)
        co = [252 159 91]/255;
        if(i == idx_SXT)
            w = 4;
            scatter(cosd(angle_SXT(i))*(SXT_reshape(i)+0.03), sind(angle_SXT(i))*(SXT_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    elseif(i >15 && i <= 30)
        co = [86 238 244]/255;
        if(i == idx_SXT+15)
            w = 4;
            scatter(cosd(angle_SXT(i))*(SXT_reshape(i)+0.03), sind(angle_SXT(i))*(SXT_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    else
        co = [50 232 117]/255;
        if(i == idx_SXT+30)
            w = 4;
            scatter(cosd(angle_SXT(i))*(SXT_reshape(i)+0.03), sind(angle_SXT(i))*(SXT_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    end
     h = plot([0 x_SXT(i)], [0 y_SXT(i)], 'linewidth', w, 'Color', co);
end

for k1 = 1:3
    plot(x_SXT_tnr(:, k1), y_SXT_tnr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_SXT_tpr(:, k1), y_SXT_tpr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_SXT_acc(:, k1), y_SXT_acc(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
end

for i = 1:length(x_GM)
    if(i <= 15)
        co = [252 159 91]/255;
        if(i == idx_GM)
            w = 4;
            scatter(cosd(angle_GM(i))*(GM_reshape(i)+0.03), sind(angle_GM(i))*(GM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    elseif(i >15 && i <= 30)
        co = [86 238 244]/255;
        if(i == idx_GM+15)
            w = 4;
            scatter(cosd(angle_GM(i))*(GM_reshape(i)+0.03), sind(angle_GM(i))*(GM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    else
        co = [50 232 117]/255;
        if(i == idx_GM+30)
            w = 4;
            scatter(cosd(angle_GM(i))*(GM_reshape(i)+0.03), sind(angle_GM(i))*(GM_reshape(i)+0.03), 30, 'r', '*');
        else
            w = 1;
        end
    end
     h = plot([0 x_GM(i)], [0 y_GM(i)], 'linewidth', w, 'Color', co);
end

for k1 = 1:3
    plot(x_GM_tnr(:, k1), y_GM_tnr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_GM_tpr(:, k1), y_GM_tpr(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
    plot(x_GM_acc(:, k1), y_GM_acc(:, k1),'-', 'linewidth', wid_arc, 'color', d(k1, :), 'LineStyle', line_t{k1})
end

x = cosd(0:30:360)';
y = sind(0:30:360)';
%intermediate lines between without arc
for i = 1:length(x)
    if(mod(i, 3) == 1)
        temp = 1.3;
        wid = 3;%thickness of axis lines
    else
        temp =  1;
        wid = 1;
    end
    plot([0 x(i)]*temp, [0 y(i)]*temp,'k-', 'linewidth', wid)
end

%line around arc - thicker
[x, y] = drawCurve(0, 360, 1.3);
plot(x,y,'-', 'linewidth', 3, 'color', [169, 217, 142]/255)
[x, y] = drawCurve(0, 270, 1.3);
plot(x,y,'-', 'linewidth', 3, 'color', [215, 241, 255]/255)
[x, y] = drawCurve(0, 180, 1.3);
plot(x,y,'-', 'linewidth', 3, 'color', [245, 211, 233]/255)
[x, y] = drawCurve(0, 90, 1.3);
plot(x,y,'-', 'linewidth', 3, 'color', [255, 217, 102]/255)
set(gca, 'fontsize', 26, 'FontWeight','bold','linewidth',4,'yticklabel',[])

