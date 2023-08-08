%% Fit YAPTAZ circadian model to stiffness data
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600];
% p0(end) = 10;
lowerLim = p0/4;
upperLim = p0*4;
lowerLim(1) = 8*3600;
upperLim(1) = 8*3600;%16*3600;
lowerLim(7) = 8*3600;%4*3600;
upperLim(7) = 8*3600;
fixedParam = [2, 8, 14, 17, 22, 25];
for i = fixedParam
    lowerLim(i) = p0(i);
    upperLim(i) = p0(i);
end
expParam = [2,8,14,17,22,25]; 
for i = expParam
    if lowerLim(i)<1
        lowerLim(i) = 1;
    end
end
options = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,'PlotFcn','pswplotbestf');%,...
%     'FunctionTolerance',1e-4,'PlotFcn','pswplotbestf','MaxStallIterations',15);%, 'MaxTime', 60*20);
pSol = particleswarm(@pToObj_CircadianClockFibroblast, length(p0), lowerLim, upperLim, options);%p0);
pToObj_CircadianClockFibroblast(pSol)
save('pSol_var4x_DynOnly_lucModel_restrictMean_1Delay.mat','pSol')

%%
pToObj_CircadianClockFibroblast(pSol)

%% Try Bayesian approach - generate dynamics with error bars
numSamples = 1e6;
tSample = 0:7*24; % say one sample per hour for one week 
periodVec = [23.8447, 24.5874, 25.5922,...
             24.7122, 25.1024, 25.1707,...
             24.2103, 24.9957, 24.9700, 24.6352,...
             24.9389, 25.4105, 25.9694,...
             23.6788, 23.2228, 22.3938, 20.4041];
periodStdVec = [0.1748, 0.2330, 0.3350,...
                0.1659, 0.2341, 0.3024,...
                0.2060, 0.1803, 0.0258, 0.3991,...
                0.1048, 0.1921, 0.2620,...
                0.2487, 0.2902, 0.4145, 0.1658];
amplVec = [1, 2.2948, 1.9392,...
           1, 1.3063, 1.1532,...
           1, 1.2874, 1.6552, 1.9425,...
           1, 1.6118, 2.3765,...
           1, 0.9531, 1.0938, 0.2969];
amplStdVec = [0.3739, 0.1094, 0.2553,...
              0.0090, 0.1081, 0.1081,...
              0.1954, 0.1494, 0.3218, 0.3908,...
              0.0235, 0.1412, 0.0588,...
              0.1250, 0.1797, 0.2344, 0.0703];
numConditions = length(periodVec);
meanDynVals = zeros(numConditions, length(tSample));
stdDynVals = zeros(numConditions, length(tSample));
figure
hold on
for i = 1:numConditions
    randVals = randn([numSamples, 2]);
    dynMat = zeros(numSamples,length(tSample));
    for j = 1:numSamples
        periodCur = periodVec(i) + periodStdVec(i)*randVals(j,1);
        amplCur = amplVec(i) + amplStdVec(i)*randVals(j,2);
        % amplCur = 1;%amplCur/amplVec(i);
        dynMat(j,:) = 0.5*amplCur*(cos(2*pi*tSample/(periodCur)));
        % plot(tSample, dynMat(j,:))
    end
    meanDynVals(i,:) = mean(dynMat,1);
    stdDynVals(i,:) = std(dynMat,1);
    % figure
    errorbar(tSample, meanDynVals(i,:), stdDynVals(i,:), stdDynVals(i,:))
    % plot(tSample, 0.5*amplVec(i)*(cos(2*pi*tSample/(periodVec(i)))))
    plot(tSample, 0.5*(cos(2*pi*tSample/(periodVec(i)))))
end
shiftVal = min(meanDynVals(:) - stdDynVals(:));
dataVec = [meanDynVals-shiftVal; stdDynVals];
y = struct;
y.dataVec = dataVec;
fixedParam = [2, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
% logLikelihood = mechanoCircadian_logLikelihood(pSol(varyLogic)./p0(varyLogic)', y)
pAll = ones(3,23);
%%
logLikelihood = mechanoCircadian_logLikelihoodDyn(pSol(varyLogic)'./p0(varyLogic)', y)

%% test with period vec
dataVec = [periodVec'; periodStdVec'];
dataVec(:,2) = [amplVec'; amplStdVec'];
y = struct;
y.dataVec = dataVec;
fixedParam = [2, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
pAll = ones(3,23);
logLikelihood = mechanoCircadian_logLikelihoodPeriod(pAll, y)

%%
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600;...
    100; 1; 10; 2];
rng(100,'twister')
uqlab
paramNames = {'tauB','nB','KeB0','KiB','KdBP','KdB',...
              'tauP','nP','KeP0','KaP','KdP',...
              'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
              'Kcap', 'StiffThresh',...
              'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
              'LatBSens','JasMag','JasSens','KdLuc',...
              'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim'};
numParam = 33;
PriorOpts.Name = 'Prior distribution on mechano-Circadian parameters';
fixedParam = [1, 2, 7, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
for i = 1:numParam
    PriorOpts.Marginals(i).Name = paramNames{i};
    if any(i==fixedParam)
        PriorOpts.Marginals(i).Type = 'Constant';
        PriorOpts.Marginals(i).Parameters = 1;
    elseif i==1 || i==7
        PriorOpts.Marginals(i).Type = 'Uniform';
        PriorOpts.Marginals(i).Parameters = [.8 1.25];
    else
        PriorOpts.Marginals(i).Type = 'Uniform';
        PriorOpts.Marginals(i).Parameters = [.25 4];
    end
end
myPriorDist = uq_createInput(PriorOpts);
myData.y.dataVec = dataVec;
myData.Name = 'PER expression oscillation';
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = 60;
Solver.MCMC.Steps = 2000;
Solver.MCMC.Visualize.Parameters = 3;%[3 8];
Solver.MCMC.Visualize.Interval = 50;
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.LogLikelihood = @mechanoCircadian_logLikelihoodDyn;
BayesOpts.Solver = Solver;
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% convert Bayesian analysis to pSol (mean values)
meanVals = myBayesianAnalysis.Results.PostProc.Percentiles.Mean;
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP')
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
pSol = p0;
pSol(varyLogic) = pSol(varyLogic) .* modeVals';
save('bayesianAnalysis_DynOnly_fixedDelays_withMRTFParams.mat','myBayesianAnalysis')

%% plot posteriors
fixedParam = [1, 2, 7, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
paramMap = find(varyLogic);
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600;...
    100; 1; 10; 2];
figure
conversionFactors = [1/3600, 1, 3600, 1, 3600, 3600, 1/3600, 1, 3600, 1, 3600,...
    3600, 1, 1, 3600, 1, 1, 1, 1, 3600, 1, 1, 3600, 1, 1, 1, 1, 1, 3600,...
    1, 1, 1, 1];
for i = 1:length(paramMap)
    subplot(3,9,i)
    paramNum = i;
    curSamples = myBayesianAnalysis.Results.Sample(end,paramNum,:);
    [pdfCur, pCur] = ksdensity(curSamples(:)*p0(paramMap(paramNum))*conversionFactors(paramMap(paramNum)),...
        'Support','positive','BoundaryCorrection','reflection');
    plot(pCur, pdfCur)
    % hold on
    % pPSOCur = pSolPSO(paramMap(paramNum))*conversionFactors(paramMap(paramNum));
    % plot([pPSOCur pPSOCur], [0, max(pdfCur)])
    xlabel(paramNames(paramMap(paramNum)))
end

%% sample from MCMC samples
samples = myBayesianAnalysis.Results.Sample(1001:end,:,:);
numParam = size(samples,2);
numSamplesNew = 400;
popParam = ones(numSamplesNew, 33);
fixedParam = [1, 2, 7, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600;...
    100; 1; 10; 2];
flatSamples = zeros(numParam, size(samples,1)*size(samples,3));
for i = 1:size(samples,1)
    flatSamples(:,size(samples,3)*(i-1)+1:size(samples,3)*i) = samples(i,:,:);
end
rVals = randi(size(flatSamples,2), numSamplesNew, 1);
for i = 1:numSamplesNew
    popParam(i,varyLogic) = flatSamples(:,rVals(i));
end
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP')
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
pSolMAP = ones(1, 33);
pSolMAP(varyLogic) = modeVals;
popParam = [popParam; pSolMAP];
% for i = 1:numParam
%     % curSamples = samples(1,i,:);
%     % [cdfCur, pCur] = ksdensity(curSamples(:),'Function','cdf',...
%     %     'Support','positive','BoundaryCorrection','reflection');
%     % cdfCur = [0, cdfCur, 1];
%     % pCur = [0, pCur, 1.01*max(pCur)];
%     % rCur = rand(numSamplesNew, 1);
%     % [cdfCur, idx] = unique(cdfCur);
%     % pCur = pCur(idx);
%     % interpVal = interp1(cdfCur, pCur, rCur);
%     % if any(isnan(interpVal)) || any(interpVal < 0)
%     %     fprintf('inspect')
%     % end
%     % popParam(:,paramMap(i)) = interpVal;
% end
popParam = popParam .* p0';

%% load stored data
load('FibroblastStiffnessPeriod.mat')
load('FibroblastStiffnessAmpl.mat')
load('FibroblastROCKInhibPeriod.mat')
load('FibroblastROCKInhibAmpl.mat')
load('FibroblastCytDPeriod.mat')
load('FibroblastCytDAmpl.mat')
load('FibroblastJasPeriod.mat')
load('FibroblastJasAmpl.mat')
load('FibroblastLatBPeriod.mat')
load('FibroblastLatBAmpl.mat')
meanFibroblastControlPeriod = (4*FibroblastStiffnessPeriod(1,1) +...
        3*FibroblastROCKInhibPeriod(1,1) + 3*FibroblastCytDPeriod(1,1) +...
        3*FibroblastLatBPeriod(1,1) + 3*FibroblastJasPeriod(1,1))/16;
stdFibroblastControlPeriod = sqrt((4*FibroblastStiffnessPeriod(1,2) +...
        3*FibroblastROCKInhibPeriod(1,2) + 3*FibroblastCytDPeriod(1,2) +...
        3*FibroblastLatBPeriod(1,2) + 3*FibroblastJasPeriod(1,2))/16);
meanFibroblastControlAmpl = (4*FibroblastStiffnessAmpl(1,1) +...
        3*FibroblastROCKInhibAmpl(1,1) + 3*FibroblastCytDAmpl(1,1) +...
        3*FibroblastLatBAmpl(1,1) + 3*FibroblastJasAmpl(1,1))/16;
stdFibroblastControlAmpl = sqrt((4*FibroblastStiffnessAmpl(1,2) +...
        3*FibroblastROCKInhibAmpl(1,2) + 3*FibroblastCytDAmpl(1,2) +...
        3*FibroblastLatBAmpl(1,2) + 3*FibroblastJasAmpl(1,2))/16);
FibroblastStiffnessPeriod(1,:) = [meanFibroblastControlPeriod, stdFibroblastControlPeriod];
FibroblastStiffnessAmpl(1,:) = [meanFibroblastControlAmpl, stdFibroblastControlAmpl];
FibroblastROCKInhibPeriod(1,:) = [meanFibroblastControlPeriod, stdFibroblastControlPeriod];
FibroblastROCKInhibAmpl(1,:) = [meanFibroblastControlAmpl, stdFibroblastControlAmpl];
FibroblastCytDPeriod(1,:) = [meanFibroblastControlPeriod, stdFibroblastControlPeriod];
FibroblastCytDAmpl(1,:) = [meanFibroblastControlAmpl, stdFibroblastControlAmpl];
FibroblastLatBPeriod(1,:) = [meanFibroblastControlPeriod, stdFibroblastControlPeriod];
FibroblastLatBAmpl(1,:) = [meanFibroblastControlAmpl, stdFibroblastControlAmpl];
FibroblastJasPeriod(1,:) = [meanFibroblastControlPeriod, stdFibroblastControlPeriod];
FibroblastJasAmpl(1,:) = [meanFibroblastControlAmpl, stdFibroblastControlAmpl];

%% Compare stiffness data
paramMat = popParam;
expDataCell = {FibroblastStiffnessPeriod, FibroblastStiffnessAmpl};
stiffnessTests = [logspace(0, 7, 20)];
testsMat = zeros(5,length(stiffnessTests));
testsMat(1,:) = stiffnessTests;
expStiffness = [1e7, 300, 19];
plotDynMat = zeros(5,length(expStiffness));
plotDynMat(1,:) = expStiffness;
legendEntries = {'Control (glass)','300 kPa PEKK scaffold','19 kPa PCL scaffold'};
[periodStiffness, amplStiffness] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

% Compare ROCK data
paramMat = popParam;
expDataCell = {FibroblastROCKInhibPeriod, FibroblastROCKInhibAmpl};
ROCKInhibTests = 0:2:30;
testsMat = zeros(5,length(ROCKInhibTests));
testsMat(1,:) = 1e7;
testsMat(2,:) = ROCKInhibTests;
FibroblastROCKInhib = [0, 10, 20];
plotDynMat = zeros(5,length(FibroblastROCKInhib));
plotDynMat(1,:) = 1e7;
plotDynMat(2,:) = FibroblastROCKInhib;
legendEntries = {'Control','10uM Y27632','20uM Y27632'};
[periodROCKInhib, amplROCKInhib] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

% CytD Tests
paramMat = popParam;
expDataCell = {FibroblastCytDPeriod, FibroblastCytDAmpl};
CytDTests = 0:.5:6;
testsMat = zeros(5,length(CytDTests));
testsMat(1,:) = 1e7;
testsMat(3,:) = CytDTests;
FibroblastCytD = [0, 1, 2, 5];
plotDynMat = zeros(5,length(FibroblastCytD));
plotDynMat(1,:) = 1e7;
plotDynMat(3,:) = FibroblastCytD;
legendEntries = {'Control','1uM CytD','2uM CytD','5uM CytD'};
[periodCytD, amplCytD] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

%% LatB Tests
paramMat = popParam;
expDataCell = {FibroblastLatBPeriod, FibroblastLatBAmpl};
LatBTests = 0:.25:2.5;
testsMat = zeros(5,length(LatBTests));
testsMat(1,:) = 1e7;
testsMat(4,:) = LatBTests;
FibroblastLatB = [0, 1, 2];
plotDynMat = zeros(5,length(FibroblastLatB));
plotDynMat(1,:) = 1e7;
plotDynMat(4,:) = FibroblastLatB;
legendEntries = {'Control','1uM LatB','2uM LatB'};
[periodLatB, amplLatB] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

% Jas Tests
paramMat = popParam;
expDataCell = {FibroblastJasPeriod, FibroblastJasAmpl};
JasTests = 0:.05:.6;
testsMat = zeros(5,length(JasTests));
testsMat(1,:) = 1e7;
testsMat(5,:) = JasTests;
FibroblastJas = [0, 0.1, 0.2, 0.5];
plotDynMat = zeros(5,length(FibroblastJas));
plotDynMat(1,:) = 1e7;
plotDynMat(5,:) = FibroblastJas;
legendEntries = {'Control','0.1uM Jas','0.2uM Jas','0.5uM Jas'};
[periodJas, amplJas] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

%% Total summary
figure
testsCell = {stiffnessTests, ROCKInhibTests, CytDTests, LatBTests, JasTests};
periodCell = {periodStiffness, periodROCKInhib, periodCytD, periodLatB, periodJas};
amplCell = {amplStiffness, amplROCKInhib, amplCytD, amplLatB, amplJas};
expCasesCell = {expStiffness, FibroblastROCKInhib, FibroblastCytD, FibroblastLatB, FibroblastJas};
expPeriodCell = {FibroblastStiffnessPeriod, FibroblastROCKInhibPeriod,...
    FibroblastCytDPeriod, FibroblastLatBPeriod, FibroblastJasPeriod};
expAmplCell = {FibroblastStiffnessAmpl, FibroblastROCKInhibAmpl,...
    FibroblastCytDAmpl, FibroblastLatBAmpl, FibroblastJasAmpl};
xlabelCell = {'Substrate stiffness (kPa)', 'Y27632 concentration (uM)',...
    'CytD concentration (uM)','LatB concentration (uM)', 'Jas concentration (uM)'};
for i = 1:length(testsCell)
    subplot(2,5,i)
    zeroVec = zeros(length(testsCell{i}),1);
    avgPeriod = zeroVec; upperBarPeriod = zeroVec; lowerBarPeriod = zeroVec;
    avgAmpl = zeroVec; upperBarAmpl = zeroVec; lowerBarAmpl = zeroVec;
    modePeriod = zeroVec; modeAmpl = zeroVec;
    if i == 1
        normAmpl = median(amplCell{i}{end});
    else
        normAmpl = median(amplCell{i}{1});
    end
    for j = 1:length(periodCell{i})
        periodPercentiles = prctile(periodCell{i}{j}(1:end-1),[25,50,75]);
        avgPeriod(j) = periodPercentiles(2);%median(periodStiffness{i});
        upperBarPeriod(j) = periodPercentiles(3) - avgPeriod(j);
        lowerBarPeriod(j) = avgPeriod(j) - periodPercentiles(1);
        modePeriod(j) = periodCell{i}{j}(end);
        amplPercentiles = prctile(amplCell{i}{j}(1:end-1)/normAmpl,[25,50,75]);
        avgAmpl(j) = amplPercentiles(2);%median(amplStiffness{i});
        upperBarAmpl(j) = amplPercentiles(3) - avgAmpl(j);
        lowerBarAmpl(j) = avgAmpl(j) - amplPercentiles(1);
        if i==1
            modeAmpl(j) = amplCell{i}{j}(end)/amplCell{i}{end}(end);
        else
            modeAmpl(j) = amplCell{i}{j}(end)/amplCell{i}{1}(end);
        end
    end
    if length(periodCell{i}{1})==1
        plot(testsCell{i}, avgPeriod/3600, 'Color', 'b', 'LineWidth', 1)
    else
        errorbar(testsCell{i}, avgPeriod/3600, lowerBarPeriod/3600,...
            upperBarPeriod/3600, 'Color', 'b', 'LineWidth', .5)
        hold on
        plot(testsCell{i}, modePeriod/3600, 'Color', [0 0 .4], 'LineWidth', 2)
    end
    if i==1
        set(gca,'XScale','log')
        xlim([min(testsCell{i}), max(testsCell{i})])
    else
        xlim([-.05*max(testsCell{i}), max(testsCell{i})])
    end
    hold on
    errorbar(expCasesCell{i}, expPeriodCell{i}(:,1),...
        expPeriodCell{i}(:,2), expPeriodCell{i}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
    xlabel(xlabelCell{i})
    ylabel('Period (hr)')
    if i==1
        legend('Estimated parameter distribution', 'Maximum likelihood parameter set', 'Experiments')
    end
    ylim([21.8 27])
    prettyGraph
    subplot(2,5,5+i);
    if length(periodCell{i}{1})==1
        plot(testsCell{i}, avgAmpl, 'Color', 'b', 'LineWidth', 1)
    else
        errorbar(testsCell{i}, avgAmpl, lowerBarAmpl,...
            upperBarAmpl, 'Color', 'b', 'LineWidth', .5)
        hold on
        plot(testsCell{i}, modeAmpl, 'Color', [0 0 .4], 'LineWidth', 2)
    end
    if i==1
        set(gca,'XScale','log')
        xlim([min(testsCell{i}), max(testsCell{i})])
    else
        xlim([-.05*max(testsCell{i}), max(testsCell{i})])
    end
    hold on
    errorbar(expCasesCell{i}, expAmplCell{i}(:,1),...
        expAmplCell{i}(:,2), expAmplCell{i}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
    xlabel(xlabelCell{i})
    ylabel('Normalized amplitude')
    % legend('Model','Experiments')
    prettyGraph
end

%% plot all dynamics
stiffnessTests = [1e7, 300, 19,...
                  1e7, 1e7, 1e7,...
                  1e7, 1e7, 1e7, 1e7,...
                  1e7, 1e7, 1e7,...
                  1e7, 1e7, 1e7, 1e7];
ROCKInhibTests = [0, 0,  0,...
                  0, 10, 20,...
                  0, 0,  0, 0,...
                  0, 0,  0,...
                  0, 0,  0, 0];
CytDTests = [0, 0,  0,...
             0, 10, 20,...
             0, 1, 2, 5,...
             0, 0, 0,...
             0, 0, 0, 0];
LatBTests = [0, 0, 0,...
             0, 0, 0,...
             0, 0, 0, 0,...
             0, 1, 2,...
             0, 0, 0, 0];
JasTests = [0, 0, 0,...
            0, 0, 0,...
            0, 0, 0, 0,...
            0, 0, 0,...
            0, 0.1, 0.2, 0.5];
periodVec = [23.8447, 24.5874, 25.5922,...
             24.7122, 25.1024, 25.1707,...
             24.2103, 24.9957, 24.9700, 24.6352,...
             24.9389, 25.4105, 25.9694,...
             23.6788, 23.2228, 22.3938, 20.4041];
amplVec = [1, 2.2948, 1.9392,...
           1, 1.3063, 1.1532,...
           1, 1.2874, 1.6552, 1.9425,...
           1, 1.6118, 2.3765,...
           1, 0.9531, 1.0938, 0.2969];
periodTest = zeros(size(periodVec));
amplTest = zeros(size(amplVec));
legendEntries = {'Control (glass)','300 kPa PEKK scaffold','19 kPa PCL scaffold',...
                 'Control', '10uM Y27632', '20uM Y27632',...
                 'Control', '1uM CytD', '2uM CytD', '5uM CytD',...
                 'Control', '1uM LatB', '2uM LatB',...
                 'Control', '0.1uM Jas', '0.2uM Jas', '0.5uM Jas'};
groups = {1:3, 4:6, 7:10, 11:13, 14:17};
for k = 1:length(groups)
    figure
    hold on
    prettyGraph
    for i = groups{k}
        ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/pSol(18)));
        kra_cur = 1 / (1 + (LatBTests(i)/pSol(26))) + (1 + pSol(27))*JasTests(i) / pSol(28);
        inhibVec = [kra_cur, ROCKLevel, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD]
        [periodCur, amplCur, ~, ~, rawOutput] = conditionToOutputs(pSol, stiffnessTests(i), inhibVec);
        t = rawOutput{1};
        y = rawOutput{2};
        periodTest(i) = periodCur(3);
        amplTest(i) = amplCur(3);
        [~,locs] = findpeaks(y(:,2));
        tShift = t(locs(2))/(24*3600);
        plot(t/(24*3600)-tShift,20*y(:,2),'LineWidth',1)
        xlim([0 5])
        ylim([0 20])
        xlabel('Time (days)')
        ylabel('PER nuclear concentration (nM)')
    end
    legend(legendEntries{groups{k}})
end

%% assess goodness of fit w.r.t. period
periodVec = [23.8447, 24.5874, 25.5922,...
             24.7122, 25.1024, 25.1707,...
             24.2103, 24.9957, 24.9700, 24.6352,...
             24.9389, 25.4105, 25.9694,...
             23.6788, 23.2228, 22.3938, 20.4041];
[~,periodTest] = pToObj_CircadianClockFibroblast(pSol);
gof = goodnessOfFit(periodTest', periodVec', 'NMSE')

function [obj, periodTest, amplTest] = pToObj_CircadianClockFibroblast(p)
%     maxTime = 3600*360;
    % order of tests: diff stiffness, ROCK inhibitor, CytoD, Latrunculin B,
    % Jas
    stiffnessTests = [1e7, 300, 19,...
                      1e7, 1e7, 1e7,...
                      1e7, 1e7, 1e7, 1e7,...
                      1e7, 1e7, 1e7,...
                      1e7, 1e7, 1e7, 1e7];
    ROCKInhibTests = [0, 0,  0,...
                      0, 10, 20,...
                      0, 0,  0, 0,...
                      0, 0,  0,...
                      0, 0,  0, 0];
    CytDTests = [0, 0,  0,...
                 0, 10, 20,...
                 0, 1, 2, 5,...
                 0, 0, 0,...
                 0, 0, 0, 0];
    LatBTests = [0, 0, 0,...
                 0, 0, 0,...
                 0, 0, 0, 0,...
                 0, 1, 2,...
                 0, 0, 0, 0];
    JasTests = [0, 0, 0,...
                0, 0, 0,...
                0, 0, 0, 0,...
                0, 0, 0,...
                0, 0.1, 0.2, 0.5];
    periodVec = [23.8447, 24.5874, 25.5922,...
                 24.7122, 25.1024, 25.1707,...
                 24.2103, 24.9957, 24.9700, 24.6352,...
                 24.9389, 25.4105, 25.9694,...
                 23.6788, 23.2228, 22.3938, 20.4041];
    amplVec = [1, 2.2948, 1.9392,...
               1, 1.3063, 1.1532,...
               1, 1.2874, 1.6552, 1.9425,...
               1, 1.6118, 2.3765,...
               1, 0.9531, 1.0938, 0.2969];
    periodTest = zeros(size(periodVec));
    amplTest = zeros(size(amplVec));

    modelCases = [1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17]; % only non-repeats in conditions
    modelStored = cell(size(periodVec));
    minVals = zeros(size(periodVec));
    for i = modelCases
        ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/p(18)));
        kra_cur = 1 / (1 + (LatBTests(i)/p(26))) + (1 + p(27))*JasTests(i) / p(28);
        inhibVec = [kra_cur, ROCKLevel, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
        [periodCur, amplCur, t, y] = conditionToOutputs(p, stiffnessTests(i), inhibVec);
        periodTest(i) = periodCur(3);
        amplTest(i) = amplCur(3);
        oscDynamics = y(:,3);%3600*p(9)./(1 + (p(10)./y(:,1)).^p(8));
        modelStored{i} = [t, oscDynamics];
        minVals(i) = min(oscDynamics);
        if i == 1
            meanP = trapz(t, y(:,2))/range(t);
            meanB = trapz(t, y(:,1))/range(t);
        end
    end
    periodTest([4,7,11,14]) = periodTest(1); % all control cases
    periodTest = periodTest/3600;
    amplTest([4,7,11,14]) = amplTest(1); % all control cases
    minVals([4,7,11,14]) = minVals(1);
    amplTestRaw = amplTest;
    amplTest = amplTest/amplTest(1); % normalize to control amplitude
    modelStored{4} = modelStored{1};
    modelStored{7} = modelStored{1};
    modelStored{11} = modelStored{1};
    modelStored{14} = modelStored{1};

    obj = 0;%50*(sum(4*((periodTest(1:3) - periodVec(1:3))./periodVec(1:3)).^2) + ...
              % sum(3*((periodTest(4:6) - periodVec(4:6))./periodVec(4:6)).^2) + ...
              % sum(3*((periodTest(7:10) - periodVec(7:10))./periodVec(7:10)).^2) + ...
              % sum(3*((periodTest(11:13) - periodVec(11:13))./periodVec(11:13)).^2) + ...
              % sum(3*((periodTest(14:17) - periodVec(14:17))./periodVec(14:17)).^2));% +...
%           .25*(sum(4*((amplTest(1:3) - amplVec(1:3))./amplVec(1:3)).^2) + ...
%              sum(3*((amplTest(4:6) - amplVec(4:6))./amplVec(4:6)).^2) + ...
%              sum(3*((amplTest(7:10) - amplVec(7:10))./amplVec(7:10)).^2) + ...
%              sum(3*((amplTest(11:13) - amplVec(11:13))./amplVec(11:13)).^2) + ...
%              sum(3*((amplTest(14:17) - amplVec(14:17))./amplVec(14:17)).^2));
    nVec = [4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
    % figure
    for i = 1:length(periodVec)
        yNormalized = (modelStored{i}(:,2) - min(minVals))/amplTestRaw(1);
        % yNormalized = (modelStored{i}(:,2) - min(modelStored{i}(:,2)))/amplTestRaw(i);
        tTest = modelStored{i}(:,1);
        PERTest = 0.5*amplVec(i)*(cos(2*pi*tTest/(periodVec(i)*3600))) + max(amplVec)/2;
        % PERTest = 0.5*(cos(2*pi*tTest/(periodVec(i)*3600))) + 0.5;
        obj = obj + nVec(i)*25e-4*sum(((yNormalized - PERTest)./amplVec(i)).^2);
        % obj = obj + nVec(i)*25e-4*sum(((yNormalized - PERTest)).^2);
        % subplot(2,9,i)
        % title(sprintf("%d", i))
        % plot(tTest, yNormalized)
        % hold on
        % plot(tTest, PERTest)
    end
    if (meanP/meanB) > 2
        obj = obj + (meanP/meanB - 2)^2;
    elseif (meanP/meanB) < 0.5
        obj = obj + (meanB/meanP - 2)^2;
    end
end

function logLikelihood = mechanoCircadian_logLikelihoodDyn(p, y)
%     maxTime = 3600*360;
    % order of tests: diff stiffness, ROCK inhibitor, CytoD, Latrunculin B,
    % Jas
    dataVec = y.dataVec;
    p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
        0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600;...
        100; 1; 10; 2];
    fixedParam = [1, 2, 7, 8, 14, 17, 22, 25];
    varyLogic = true(length(p0),1);
    varyLogic(fixedParam) = false;
    logLikelihood = zeros(size(p,1),1);
    parfor k = 1:size(p,1)
        meanB = [];
        meanP = [];
        pCur = p0;
        pCur(varyLogic) = p(k,:).*pCur(varyLogic)';
        % pCur(7) = pCur(1);
        stiffnessTests = [1e7, 300, 19,...
            1e7, 1e7, 1e7,...
            1e7, 1e7, 1e7, 1e7,...
            1e7, 1e7, 1e7,...
            1e7, 1e7, 1e7, 1e7];
        ROCKInhibTests = [0, 0,  0,...
            0, 10, 20,...
            0, 0,  0, 0,...
            0, 0,  0,...
            0, 0,  0, 0];
        CytDTests = [0, 0,  0,...
            0, 10, 20,...
            0, 1, 2, 5,...
            0, 0, 0,...
            0, 0, 0, 0];
        LatBTests = [0, 0, 0,...
            0, 0, 0,...
            0, 0, 0, 0,...
            0, 1, 2,...
            0, 0, 0, 0];
        JasTests = [0, 0, 0,...
            0, 0, 0,...
            0, 0, 0, 0,...
            0, 0, 0,...
            0, 0.1, 0.2, 0.5];
        expMeanDyn = dataVec(1:length(JasTests), :);
        expStdDyn = dataVec(length(JasTests)+1:end, :);
        amplTest = zeros(size(JasTests));

        modelCases = [1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17]; % only non-repeats in conditions
        modelStored = cell(size(JasTests));
        minVals = zeros(size(JasTests));
        for i = modelCases
            ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/pCur(18)));
            kra_cur = 1 / (1 + (LatBTests(i)/pCur(26))) + (1 + pCur(27))*JasTests(i) / pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD]
            [~, amplCur, t, y] = conditionToOutputs(pCur, stiffnessTests(i), inhibVec);
            amplTest(i) = amplCur(3);
            oscDynamics = y(:,3);
            modelStored{i} = [t, oscDynamics];
            minVals(i) = min(oscDynamics);
            if i == 1
                meanP = trapz(t, y(:,2))/range(t);
                meanB = trapz(t, y(:,1))/range(t);
            end
        end
        amplTest([4,7,11,14]) = amplTest(1); % all control cases
        amplTestRaw = amplTest;
        % amplTest = amplTest/amplTest(1); % normalize to control amplitude
        modelStored{4} = modelStored{1};
        modelStored{7} = modelStored{1};
        modelStored{11} = modelStored{1};
        modelStored{14} = modelStored{1};

        logLikelihood(k) = 0;
        % nVec = [4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
        % figure
        for i = 1:length(JasTests)
            yNormalized = (modelStored{i}(:,2) - min(minVals))/amplTestRaw(1);
            % yNormalized = (modelStored{i}(:,2) - min(modelStored{i}(:,2)))/amplTestRaw(i);
            tTest = modelStored{i}(:,1);
            curMeanDyn = expMeanDyn(i,1:length(tTest))';
            curStdDyn = expStdDyn(i,1:length(tTest))';
            useLogic = curStdDyn > 0; % in certain cases, this may be zero
            logLikelihood(k) = logLikelihood(k) - sum(((yNormalized(useLogic) - curMeanDyn(useLogic))./curStdDyn(useLogic)).^2);
            % subplot(2,9,i)
            % title(sprintf("%d", i))
            % plot(tTest/(24*3600), yNormalized)
            % hold on
            % errorbar(tTest/(24*3600), curMeanDyn, curStdDyn, curStdDyn)
        end
        % if (meanP/meanB) > 2
        %     logLikelihood(k) = logLikelihood(k) - (meanP/meanB - 2)^2;
        % elseif (meanP/meanB) < 0.5
        %     logLikelihood(k) = logLikelihood(k) - (meanB/meanP - 2)^2;
        % end
    end
end

function logLikelihood = mechanoCircadian_logLikelihoodPeriod(p, y)
%     maxTime = 3600*360;
    % order of tests: diff stiffness, ROCK inhibitor, CytoD, Latrunculin B,
    % Jas
    dataVec = y.dataVec;
    p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; 0.4/3600;...
        0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600];
    fixedParam = [2, 8, 14, 17, 22, 25];
    varyLogic = true(length(p0),1);
    varyLogic(fixedParam) = false;
    logLikelihood = zeros(size(p,1),1);
    stiffnessTests = [1e7, 300, 19,...
        1e7, 1e7, 1e7,...
        1e7, 1e7, 1e7, 1e7,...
        1e7, 1e7, 1e7,...
        1e7, 1e7, 1e7, 1e7];
    ROCKInhibTests = [0, 0,  0,...
        0, 10, 20,...
        0, 0,  0, 0,...
        0, 0,  0,...
        0, 0,  0, 0];
    CytDTests = [0, 0,  0,...
        0, 10, 20,...
        0, 1, 2, 5,...
        0, 0, 0,...
        0, 0, 0, 0];
    LatBTests = [0, 0, 0,...
        0, 0, 0,...
        0, 0, 0, 0,...
        0, 1, 2,...
        0, 0, 0, 0];
    JasTests = [0, 0, 0,...
        0, 0, 0,...
        0, 0, 0, 0,...
        0, 0, 0,...
        0, 0.1, 0.2, 0.5];
    expMeanPeriod = dataVec(1:length(JasTests), 1);
    expStdPeriod = dataVec(length(JasTests)+1:end, 1);
    expMeanAmpl = dataVec(1:length(JasTests), 2);
    expStdAmpl = dataVec(length(JasTests)+1:end, 2);
    for k = 1:size(p,1)
        meanP = [];
        meanB = [];
        pCur = p0;
        pCur(varyLogic) = p(k,:).*pCur(varyLogic)';
        amplTest = zeros(size(JasTests));
        periodTest = zeros(size(JasTests));
        modelCases = [1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17]; % only non-repeats in conditions
        modelStored = cell(size(JasTests));
        minVals = zeros(size(JasTests));
        for i = modelCases
            ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/pCur(18)));
            kra_cur = 1 / (1 + (LatBTests(i)/pCur(26))) + (1 + pCur(27))*JasTests(i) / pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD]
            [periodCur, amplCur, t, y] = conditionToOutputs(pCur, stiffnessTests(i), inhibVec);
            amplTest(i) = amplCur(3);
            periodTest(i) = periodCur(3);
            oscDynamics = y(:,3);
            modelStored{i} = [t, oscDynamics];
            minVals(i) = min(oscDynamics);
            if i == 1
                meanP = trapz(t, y(:,2))/range(t);
                meanB = trapz(t, y(:,1))/range(t);
            end
        end
        amplTest([4,7,11,14]) = amplTest(1); % all control cases
        amplTest = amplTest/amplTest(1);
        periodTest([4,7,11,14]) = periodTest(1); % all control cases
        periodTest = periodTest/3600;   
        logLikelihood(k) = -sum(((periodTest'-expMeanPeriod)./expStdPeriod).^2) -...
            0.1*sum(((amplTest'-expMeanAmpl)./expStdAmpl).^2);
        if (meanP/meanB) > 2
            logLikelihood(k) = logLikelihood(k) - (meanP/meanB - 2)^2;
        elseif (meanP/meanB) < 0.5
            logLikelihood(k) = logLikelihood(k) - (meanB/meanP - 2)^2;
        end
    end
end

function [periodTests, amplTests] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries)
% testsMat has numTests columns and 6 rows: stiffness, ROCK inhib, CytD, LatB, Jas
% plotDynMat has similar structure, but contains entries to plot (number of
% curves given by number of columns in the matrix)
    
    periodTests = cell(size(testsMat,2),1);
    amplTests = cell(size(testsMat,2),1);
    for k = 1:size(paramMat,1)
        pCur= paramMat(k,:);
        for i = 1:size(testsMat,2)
            ROCKLevel =  1 / (1 + (testsMat(2,i)/pCur(18)));
            kra_cur = 1 / (1 + (testsMat(4,i)/pCur(26))) +...
                (1 + pCur(27))*testsMat(5,i)/pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 1, testsMat(3,i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
            [curPeriod, curAmpl] = conditionToOutputs(pCur, testsMat(1,i), inhibVec, 3600*480);
            periodTests{i}(k) = curPeriod(3);
            amplTests{i}(k) = curAmpl(3);
        end
    end

    figure
    sp3 = subplot(1,3,3);
    xlabel('Time (days)')
    ylabel('PER:Luc luminescence (au)')
    prettyGraph
    hold on
    oscDynamics = cell(size(plotDynMat,2));
    % figure
    for k = 1:size(paramMat,1)
        pCur= paramMat(k,:);
        for i = 1:size(plotDynMat,2)
            ROCKLevel =  1 / (1 + (plotDynMat(2,i)/pCur(18)));
            kra_cur = 1 / (1 + (plotDynMat(4,i)/pCur(26))) +...
                (1+pCur(27))*plotDynMat(5,i) / pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 1, plotDynMat(3,i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
            [~, ~, t, y] = conditionToOutputs(pCur, plotDynMat(1,i), inhibVec, 3600*480);
            storeLogic = t <= 24*5*3600;
            oscDynamics{i}(k,:) = y(storeLogic,3);
            % plot(t,y(:,3))
            % hold on
            % drawnow
        end
    end
    tStored = t(storeLogic);
    meanDyn = zeros(size(plotDynMat,2),length(tStored));
    stdDyn = zeros(size(plotDynMat,2),length(tStored));
    popLogic = size(paramMat,1) > 1;
    for i = 1:size(plotDynMat,2)
        meanDyn(i,:) = mean(oscDynamics{i}(1:end-1,:),1);
        if ~popLogic
            plot(tStored/(24*3600), meanDyn(i,:),'LineWidth',1)
        else
            stdDyn(i,:) = std(oscDynamics{i}(1:end-1,:),1)/sqrt(size(oscDynamics{i},1));
            errorbar(tStored/(24*3600), meanDyn(i,:), stdDyn(i,:), stdDyn(i,:))
        end
    end
    xlim([0 5])
    legend(legendEntries)

    avgPeriod = zeros(size(testsMat,2),1);
    upperBarPeriod = zeros(size(testsMat,2),1);
    lowerBarPeriod = zeros(size(testsMat,2),1);
    avgAmpl = zeros(size(testsMat,2),1);
    upperBarAmpl = zeros(size(testsMat,2),1);
    lowerBarAmpl = zeros(size(testsMat,2),1);
    for i = 1:size(testsMat,2)
        periodPercentiles = prctile(periodTests{i}(1:end-1),[25,50,75]);
        avgPeriod(i) = periodTests{i}(end);%median(periodStiffness{i});
        % stdPeriod(i) = std(periodStiffness{i})/sqrt(length(periodStiffness{i}));
        upperBarPeriod(i) = periodPercentiles(3) - avgPeriod(i);
        lowerBarPeriod(i) = avgPeriod(i) - periodPercentiles(1);
        amplPercentiles = prctile(amplTests{i}(1:end-1),[25,50,75]);
        avgAmpl(i) = amplTests{i}(end);%median(amplStiffness{i});
        % stdAmpl(i) = std(amplStiffness{i})/sqrt(length(periodStiffness{i}));
        upperBarAmpl(i) = amplPercentiles(3) - avgAmpl(i);
        lowerBarAmpl(i) = avgAmpl(i) - amplPercentiles(1);
    end
    sp1 = subplot(1,3,1);
    activeTestsLogic = false(size(testsMat,1),1);
    for i = 1:size(testsMat,1)
        activeTestsLogic(i) = ~all(testsMat(i,:)==testsMat(i,1)); 
    end
    testStrings = {'Substrate stiffness (kPa)', 'Y27632 concentration (uM)'...
        'CytD concentration (uM)', 'LatB concentration (uM)', 'Jasplakinolide concentration (uM)'};
    xlog = false;
    controlIdx = 1;
    if sum(activeTestsLogic)==1
        testsVec = testsMat(activeTestsLogic,:);
        expVec = plotDynMat(activeTestsLogic,:);
        xlabelString = testStrings{activeTestsLogic};
        if find(activeTestsLogic)==1
            xlog = true;
            controlIdx = size(testsMat,2);
        end
    else
        testsVec = 1:size(testsMat,2);
        expVec = 1:size(plotDynMat,2);
        xlabelString = 'Test number';
    end
    if ~popLogic
        plot(testsVec, avgPeriod/3600, 'Color', 'b', 'LineWidth', 1)
    else
        errorbar(testsVec, avgPeriod/3600, lowerBarPeriod/3600, upperBarPeriod/3600,...
            'Color', 'b', 'LineWidth', 1)
    end
    if xlog
        set(gca, 'XScale', 'log')
    end
    hold on
    errorbar(expVec, expDataCell{1}(:,1),...
        expDataCell{1}(:,2), expDataCell{1}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
    xlabel(xlabelString)
    ylabel('Period (hr)')
    legend('Model','Experiments')
    ylim([20 28])
    prettyGraph
    sp2 = subplot(1,3,2);
    if ~popLogic
        plot(testsVec, avgAmpl/avgAmpl(controlIdx), 'Color', 'b', 'LineWidth', 1)
    else
        errorbar(testsVec, avgAmpl/avgAmpl(controlIdx), lowerBarAmpl/avgAmpl(controlIdx),...
            upperBarAmpl/avgAmpl(controlIdx), 'Color', 'b', 'LineWidth', 1)
    end
    if xlog
        set(gca, 'XScale', 'log')
    end
    hold on
    errorbar(expVec, expDataCell{2}(:,1),...
        expDataCell{2}(:,2), expDataCell{2}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
    xlabel(xlabelString)
    ylabel('Normalized amplitude')
    legend('Model','Experiments')
    ylim([0 3])
    prettyGraph
    pos1 = get(sp1, 'Position');
    pos2 = get(sp2, 'Position');
    pos3 = get(sp3, 'Position');
    set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
    set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
    set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])
end