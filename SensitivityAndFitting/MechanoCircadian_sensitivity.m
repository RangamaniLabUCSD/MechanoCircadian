%% Driver script for MechanoCircadian model sensitivity analysis using UQLab
clearvars
rng(100,'twister')
uqlab
ModelOpts.mFile = 'uq_MechanoCircadian';
myModel = uq_createModel(ModelOpts);
paramNames = {'tauB','nB','KeB0','KiB','KdBP','KdB',...
              'tauP','nP','KeP0','KaP','KdP',...
              'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
              'ROCKInhibSens', 'StiffThresh',...
              'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
              'LatBSens','JasMag','JasSens','KdLuc',...
              'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim','Kcap'};
activeIdx = [1,3:7,9:13, 15,16,19,20,21,23,24,29,30,31,32]; % exclude Hill coefficients and parameters associated with inhibitor treatments
paramNames = paramNames(activeIdx);
numParam = 22;
for i = 1:numParam
    InputOpts.Marginals(i).Type = 'Lognormal';
    InputOpts.Marginals(i).Name = paramNames{i};
    if i==1 || i==7
        InputOpts.Marginals(i).Parameters = [0 log(2)];
    else
        InputOpts.Marginals(i).Parameters = [0 log(4)];
    end
end
myInput = uq_createInput(InputOpts);
SobolOpts.Type = 'uq_sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 1000;
SobolOpts.Bootstrap.Replications = 100;
SobolOpts.Bootstrap.Alpha = 0.05;
SobolAnalysis = uq_createAnalysis(SobolOpts);
save('SobolAnalysisNew_1000_newLucDecay.mat','SobolAnalysis')
% display
uq_visualizeSobolIndices(SobolAnalysis.Results, [3, 6])

%% customized plot for Figure 3
plotIdx = 6;
curColor = 'r';
figure
% subplot(1,2,1)
hold on
ax = gca;
uq_formatDefaultAxes(ax)
% Create the bar plot
Results = SobolAnalysis.Results;
keepIndices = [1:13,15:size(Results.Total,1)];
NumOfIndices = length(keepIndices);
x = 1:NumOfIndices;
curData = Results.Total(keepIndices,plotIdx);
[curData, sortIdx] = sort(curData, 'descend');
uq_bar(x, curData,'FaceColor', curColor);
% Add error bars if Bootstrap option is set
if isfield(Results, 'Bootstrap')
    meanData = Results.Bootstrap.Total.Mean(keepIndices,plotIdx);
    lb = meanData - Results.Bootstrap.Total.CI(keepIndices,plotIdx,1);
    ub = Results.Bootstrap.Total.CI(keepIndices,plotIdx,2) - meanData;
    meanData = meanData(sortIdx);
    lb = lb(sortIdx);
    ub = ub(sortIdx);
    errorbar(x, meanData, lb, ub,...
        'kx', 'LineStyle', 'none', 'LineWidth', 1);
end
xticks(x)
paramNamesCur = paramNames(keepIndices);
paramNamesCur = paramNamesCur(sortIdx);
xticklabels(paramNamesCur)
ylabel('Total Sobol indices')
prettyGraph