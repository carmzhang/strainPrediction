%%
% Resistance prediction using whole genome sequencing through three
% approaches - a literature search, CARD database, and resFinder database.
% Choose the drug and database of interest in line 12.
% Supplementary Table 3.2
%
%%

clear
close all

%% CHOOSE drug and database for resistance prediction
% drug - OPTIONS: 1 - SAM
%                 2 - GM
%                 3 - SXT
%                 4 - CIP
% database - OPTIONS: 'lit' - literature search
%                     'CARD' - CARD database
%                     'resfinder' - resFinder database
% EX: chosen drug is SAM and databse is literature search
drug = 1;
database = 'lit';

if(strcmp(database, 'lit'))
    [sens, spec, acc] = resistance_literature(drug);
elseif(strcmp(database, 'card'))
    [sens, spec, acc] = resistance_card(drug);
elseif(strcmp(database, 'resfinder'))
    [sens, spec, acc] = resistance_resfinder(drug);
end

disp('Sensitivity: ')
disp(sens)
disp('Specificity: ')
disp(spec)
disp('Accuracy: ')
disp(acc)
