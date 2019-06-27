function [sens, spec, acc] = resistance_literature(drug)
%FUNCTION: returns the sensitivity, specificity, and accuracy of the 
%          resistance prediction using the resFinder database
% INPUT:
%   drug - which drug to report resFinder prediction for
%           OPTIONS: 1 - SAM
%                    2 - GM
%                    3 - SXT
%                    4 - CIP
%
% OUTPUT: metrics for resistance prediction
%   sens - sensitivity
%   spec - specificity
%   acc - accuracy

%%
% Resistance prediction using literature search of genetic basis for
% resistance.
%
%%

%actual phenotype: SAM, GM, SXT, CIP

load('resistance_WGSAnalysis_20180501.mat') %with 95% similarity and 50% gene length

%remove whitespace in names
sNames_p = strrep(sNames_phenotypic, ' ', '');
sNames_w = strrep(sNames_WGS, ' ', '');
setdiff(sNames_p, sNames_w)

threshold = 95; %for acove only

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

comb = drugProfile + predictedProfile;%SAM, GM, SXT, CIP
minus = drugProfile - predictedProfile;

oneTP = sum(comb(:, drug) == 2);
oneTN = sum(comb(:, drug) == 0);
oneFP = sum(minus(:, drug) == -1);
oneFN = sum(minus(:, drug) == 1);

sens = oneTP/sum(drugProfile(:, drug) == 1);
spec = oneTN/sum(drugProfile(:, drug) == 0);
acc = (oneTP+oneTN)/length(sNames_phenotypic);




end

