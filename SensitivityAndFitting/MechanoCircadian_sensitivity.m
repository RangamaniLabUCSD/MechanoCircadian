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
    InputOpts.Marginals(i).Type = 'Uniform';
    InputOpts.Marginals(i).Name = paramNames{i};
    if i==1 || i==6
        InputOpts.Marginals(i).Parameters = [.8 1.2];
    else
        InputOpts.Marginals(i).Parameters = [.25 4];
    end
end
myInput = uq_createInput(InputOpts);
SobolOpts.Type = 'uq_sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 10000;
SobolOpts.Bootstrap.Replications = 100;
SobolOpts.Bootstrap.Alpha = 0.05;
SobolAnalysis = uq_createAnalysis(SobolOpts);
save('SobolAnalysisNew_50000_newLucDecay.mat','SobolAnalysis')
% display
uq_visualizeSobolIndices(SobolAnalysis.Results, [3, 6])

%% customized plot for Figure S1
categoryIdx = {[1:11, 29], [12:17, 20:25], [18:19, 26:28, 33:34], 30:32}; % sorting into four categories of parameters
categoryNames = {'Circadian', 'Mechano-circadian', 'Treatments', 'MRTF'};
activeIdx = [1,3:7,9:13, 15,16,19,20,21,23,24,29,30,31,32]; % exclude Hill coefficients and parameters associated with inhibitor treatments
categoryIdxSel = cell(1,4);
for i = 1:length(categoryIdx)
    curCategory = [];
    for j = 1:length(categoryIdx{i})
        curIdx = find(activeIdx == categoryIdx{i}(j));
        if ~isempty(curIdx)
            curCategory = [curCategory, curIdx]; %#ok<AGROW>
        end
    end
    categoryIdxSel{i} = curCategory;
end

plotIdx = 3; % number for qoi you want to plot
curColor = 'b';
figure
% subplot(1,2,1)
hold on
ax = gca;
uq_formatDefaultAxes(ax)
% Create bar plot sorted from highest total order Sobol index to lowest,
% sorted into categories
Results = SobolAnalysis.Results;
xCur = 1;
xAll = [];
storeNamesAll = {};
for i = [1,2,4] % not plotting any treatment parameters
    curIndices = categoryIdxSel{i};
    % keepIndices = [1:13,15:size(Results.Total,1)]; % parameters to plot
    NumOfIndices = length(curIndices);
    x = xCur:xCur+NumOfIndices-1;
    curData = Results.Total(curIndices,plotIdx);
    [curData, sortIdx] = sort(curData, 'descend');
    uq_bar(x, curData,'FaceColor', curColor);
    % Add error bars if Bootstrap option is set
    if isfield(Results, 'Bootstrap')
        meanData = Results.Bootstrap.Total.Mean(curIndices,plotIdx);
        lb = meanData - Results.Bootstrap.Total.CI(curIndices,plotIdx,1);
        ub = Results.Bootstrap.Total.CI(curIndices,plotIdx,2) - meanData;
        meanData = meanData(sortIdx);
        lb = lb(sortIdx);
        ub = ub(sortIdx);
        errorbar(x, meanData, lb, ub,...
            'kx', 'LineStyle', 'none', 'LineWidth', 1);
    end
    storeNames = paramNames(curIndices);
    storeNames = storeNames(sortIdx);
    xAll = [xAll, x]; %#ok<AGROW>
    storeNamesAll(length(storeNamesAll)+1:length(storeNamesAll)+length(x)) = storeNames;
    xCur = xCur+NumOfIndices+2;
    hold on
end
xticks(xAll)
xticklabels(storeNamesAll)
ylabel('Total Sobol indices')
prettyGraph