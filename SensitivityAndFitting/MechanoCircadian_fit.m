%% Script used for Bayesian parameter estimation of MechanoCircadian model
% Emmet Francis, 2023

%% load stored data from Xiong et al
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

%% Generate dynamics with error bars for Bayesian fitting
numSamples = 1e6;
tSample = 0:7*24; % say one sample per hour for one week 
periodVec = [FibroblastStiffnessPeriod(:,1)', FibroblastROCKInhibPeriod(:,1)',...
             FibroblastCytDPeriod(:,1)', FibroblastLatBPeriod(:,1)', FibroblastJasPeriod(:,1)'];
periodStdVec = [FibroblastStiffnessPeriod(:,2)', FibroblastROCKInhibPeriod(:,2)',...
             FibroblastCytDPeriod(:,2)', FibroblastLatBPeriod(:,2)', FibroblastJasPeriod(:,2)'];
amplVec = [FibroblastStiffnessAmpl(:,1)', FibroblastROCKInhibAmpl(:,1)',...
             FibroblastCytDAmpl(:,1)', FibroblastLatBAmpl(:,1)', FibroblastJasAmpl(:,1)'];
amplStdVec = [FibroblastStiffnessAmpl(:,2)', FibroblastROCKInhibAmpl(:,2)',...
             FibroblastCytDAmpl(:,2)', FibroblastLatBAmpl(:,2)', FibroblastJasAmpl(:,2)'];
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
end
shiftVal = min(meanDynVals(:));% - stdDynVals(:));
dataVec = [meanDynVals; stdDynVals];
y = struct;
y.dataVec = dataVec;
y.periodVec = periodVec;
y.periodStdVec = periodStdVec;

%% Bayesian parameter estimation using UQLab
p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
rng(100,'twister')
uqlab
paramNames = {'tauB','nB','KeB0','KiB','KdBP','KdB',...
              'tauP','nP0','KeP0','KaP','KdP',...
              'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
              'ROCKInhibSens', 'StiffThresh',...
              'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
              'LatBSens','JasMag','JasSens','KdLuc',...
              'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim','Kcap',...
              'KeP1', 'KiP', 'KdR', 'tauR',...
              'KeR2,Y', 'KYR', 'KeR2,M', 'KMR', 'nP1', 'nYR', 'nMR'};
numParam = length(p0);
PriorOpts.Name = 'Prior distribution on mechano-Circadian parameters';
fixedParam = [5, 14, 17, 22, 25, 44:45];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
for i = 1:numParam
    PriorOpts.Marginals(i).Name = paramNames{i};
    if any(i==fixedParam)
        PriorOpts.Marginals(i).Type = 'Constant';
        PriorOpts.Marginals(i).Parameters = 1;
    elseif i==1 || i==7 || i==38 || i==2 || i==8 || i==43
        PriorOpts.Marginals(i).Type = 'Uniform';
        PriorOpts.Marginals(i).Parameters = [.5 2];
    else
        PriorOpts.Marginals(i).Type = 'Uniform';
        PriorOpts.Marginals(i).Parameters = [.25 4];
    end
end
myPriorDist = uq_createInput(PriorOpts);
myData.y = y;%.dataVec = dataVec;
myData.Name = 'PER expression oscillation';
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = 60;
Solver.MCMC.Steps = 5000;
Solver.MCMC.Visualize.Parameters = 3;%[3 8];
Solver.MCMC.Visualize.Interval = 50;
Solver.MCMC.Seed = initialSeed;
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.LogLikelihood = @mechanoCircadian_logLikelihoodDyn;
BayesOpts.Solver = Solver;
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% convert Bayesian analysis to pSol (MAP values)
meanVals = myBayesianAnalysis.Results.PostProc.Percentiles.Mean;
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP')
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
pSol = p0;
pSol(varyLogic) = pSol(varyLogic) .* modeVals';
save('bayesianAnalysis_RBPFinal_5001to10000_varypExp_weightPeriod100.mat','myBayesianAnalysis')

%% Test convergence of MCMC 
% (code from Nathaniel Linden - https://github.com/RangamaniLabUCSD/CIUKF-MCMC/blob/main/utils/computeIACT.m)
samples = myBayesianAnalysis.Results.Sample;
numParam = size(samples,2);
numWalkers = size(samples,3);
Nsteps = size(samples,1);
IACTs = zeros(numParam,numWalkers); % matrix of IACTs (nParam x nWalkers)
for k=1:numWalkers
    for j=1:numParam
        % parameters here are essentially all default
        [~,~,~,tauinttmp,~,~] = UWerr_fft(squeeze(samples(:,j,k)),1.5,Nsteps,1,1,1);
        IACTs(j,k) = tauinttmp;
    end
end
% Take the average over each walker:
meanIACTs = mean(IACTs,2);
% Take the max over all set of params
IACT = max(meanIACTs);

%%
pTest = exp(0.5*randn([500,length(modeVals)]));
LLVals = zeros(size(pTest,1),1);
for i = 1:size(pTest,1)
    pCur = p0;
    pCur(varyLogic) = pCur(varyLogic).*pTest(i,:)';
    LLVals(i) = computeLogLikelihood(pCur, y.dataVec);
    fprintf("Computed LL %d of %d\n", i, size(pTest,1))
end

%%
% uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP')
% modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
pSol = p0;
pSol(varyLogic) = pSol(varyLogic) .* modeVals';
% [~,maxIdx] = max(LLVals);
% pSol(varyLogic) = pSol(varyLogic).*pTest(maxIdx,:)';
computeLogLikelihood(pSol, y, true)

%% plot posteriors
fixedParam = [5, 14, 17, 22, 25, 44:45];
p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
paramMap = find(varyLogic);
figure
conversionFactors = [1/3600, 1, 3600, 1, 3600, 3600, 1/3600, 1, 3600, 1, 3600,...
    3600, 1, 1, 3600, 1, 1, 1, 1, 3600, 1, 1, 3600, 1, 1, 1, 1, 1, 3600,...
    1, 1, 1, 1, 1, 3600, 1, 3600, 1/3600, 3600, 1, 3600, 1, 1, 1, 1];
for i = 1:length(paramMap)
    subplot(5,8,i)
    paramNum = i;
    curSamples = myBayesianAnalysis.Results.Sample(501:end,paramNum,:);
    curMagn = p0(paramMap(paramNum))*conversionFactors(paramMap(paramNum));
    if paramMap(i)==1 || paramMap(i)==7
        [pdfCur, pCur] = ksdensity(curSamples(:)*curMagn,...
            'Support','positive','BoundaryCorrection','reflection','Bandwidth',curMagn*.06);
    else
        [pdfCur, pCur] = ksdensity(curSamples(:)*curMagn,...
            'Support','positive','BoundaryCorrection','reflection','Bandwidth',curMagn*.2);
    end
    plot(pCur, pdfCur*curMagn,'LineWidth',1)
    hold on
    xlabel(paramNames(paramMap(paramNum)))
    prettyGraph
end

%% sample from MCMC samples
samples = myBayesianAnalysis.Results.Sample(501:end,:,:); % discard first samples as burnin
p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
numParam = size(samples,2);
numSamplesNew = 50;
popParam = ones(numSamplesNew, length(p0));
fixedParam = [5, 14, 17, 22, 25, 44:45];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
flatSamples = zeros(numParam, size(samples,1)*size(samples,3));
for i = 1:size(samples,1)
    flatSamples(:,size(samples,3)*(i-1)+1:size(samples,3)*i) = samples(i,:,:);
end
% flatSamples = exp(0.05*randn([numParam,200]));
% for i = 1:size(flatSamples,2)
%     flatSamples(:,i) = pTest(maxIdx,:)'.*flatSamples(:,i);
% end

rVals = randi(size(flatSamples,2), numSamplesNew, 1);
for i = 1:numSamplesNew
    popParam(i,varyLogic) = flatSamples(:,rVals(i));
    % popParam(i,31:32) = 10*popParam(i,31:32);
    % popParam(i,[16,24,42]) = 10*popParam(i,[16,24,42]);
end
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP')
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
pSolMAP = ones(1, length(p0));
pSolMAP(varyLogic) = modeVals;
popParam = [popParam; pSolMAP];
popParam = popParam .* p0';

%% Average out control data for plotting
meanFibroblastControlPeriod = (4*FibroblastStiffnessPeriod(1,1) +...
        3*FibroblastROCKInhibPeriod(1,1) + 3*FibroblastCytDPeriod(1,1) +...
        3*FibroblastLatBPeriod(1,1) + 3*FibroblastJasPeriod(1,1))/16;
stdFibroblastControlPeriod = sqrt((4*FibroblastStiffnessPeriod(1,2)^2 +...
        3*FibroblastROCKInhibPeriod(1,2)^2 + 3*FibroblastCytDPeriod(1,2)^2 +...
        3*FibroblastLatBPeriod(1,2)^2 + 3*FibroblastJasPeriod(1,2)^2)/16);
meanFibroblastControlAmpl = (4*FibroblastStiffnessAmpl(1,1) +...
        3*FibroblastROCKInhibAmpl(1,1) + 3*FibroblastCytDAmpl(1,1) +...
        3*FibroblastLatBAmpl(1,1) + 3*FibroblastJasAmpl(1,1))/16;
stdFibroblastControlAmpl = sqrt((4*FibroblastStiffnessAmpl(1,2)^2 +...
        3*FibroblastROCKInhibAmpl(1,2)^2 + 3*FibroblastCytDAmpl(1,2)^2 +...
        3*FibroblastLatBAmpl(1,2)^2 + 3*FibroblastJasAmpl(1,2)^2)/16);
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
stiffnessTests = [logspace(0, 7, 30)];
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
ROCKInhibTests = 0:30;
testsMat = zeros(5,length(ROCKInhibTests));
testsMat(1,:) = 1e7;
testsMat(2,:) = ROCKInhibTests;
FibroblastROCKInhib = [0, 10, 20];
plotDynMat = zeros(5,length(FibroblastROCKInhib));
plotDynMat(1,:) = 1e7;
plotDynMat(2,:) = FibroblastROCKInhib;
legendEntries = {'Control','10μM Y27632','20μM Y27632'};
[periodROCKInhib, amplROCKInhib] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

% CytD Tests
paramMat = popParam;
expDataCell = {FibroblastCytDPeriod, FibroblastCytDAmpl};
CytDTests = 0:.25:6;
testsMat = zeros(5,length(CytDTests));
testsMat(1,:) = 1e7;
testsMat(3,:) = CytDTests;
FibroblastCytD = [0, 1, 2, 5];
plotDynMat = zeros(5,length(FibroblastCytD));
plotDynMat(1,:) = 1e7;
plotDynMat(3,:) = FibroblastCytD;
legendEntries = {'Control','1μM CytD','2μM CytD','5μM CytD'};
[periodCytD, amplCytD] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

% LatB Tests
paramMat = popParam;
expDataCell = {FibroblastLatBPeriod, FibroblastLatBAmpl};
LatBTests = 0:.1:2.5;
testsMat = zeros(5,length(LatBTests));
testsMat(1,:) = 1e7;
testsMat(4,:) = LatBTests;
FibroblastLatB = [0, 1, 2];
plotDynMat = zeros(5,length(FibroblastLatB));
plotDynMat(1,:) = 1e7;
plotDynMat(4,:) = FibroblastLatB;
legendEntries = {'Control','1μM LatB','2μM LatB'};
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
legendEntries = {'Control','0.1μM Jas','0.2μM Jas','0.5μM Jas'};
[periodJas, amplJas] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries);

%% Total summary
spaghetti = true; % optional - plot individual sample curves (spaghetti plot)
figure
testsCell = {stiffnessTests, ROCKInhibTests, CytDTests, LatBTests, JasTests};
periodCell = {periodStiffness, periodROCKInhib, periodCytD, periodLatB, periodJas};
amplCell = {amplStiffness, amplROCKInhib, amplCytD, amplLatB, amplJas};
expCasesCell = {expStiffness, FibroblastROCKInhib, FibroblastCytD, FibroblastLatB, FibroblastJas};
expPeriodCell = {FibroblastStiffnessPeriod, FibroblastROCKInhibPeriod,...
    FibroblastCytDPeriod, FibroblastLatBPeriod, FibroblastJasPeriod};
expAmplCell = {FibroblastStiffnessAmpl, FibroblastROCKInhibAmpl,...
    FibroblastCytDAmpl, FibroblastLatBAmpl, FibroblastJasAmpl};
xlabelCell = {'Substrate stiffness (kPa)', 'Y27632 concentration (μM)',...
    'CytD concentration (μM)','LatB concentration (μM)', 'Jas concentration (μM)'};
for i = 1:length(testsCell)
    subplot(2,5,i)
    zeroVec = zeros(length(testsCell{i}),1);
    if i == 1
        amplVals = amplCell{i}(end,1:end-1);
    else
        amplVals = amplCell{i}(1,1:end-1);
    end
    normAmpl = median(amplVals(~isnan(amplVals)));
    modePeriod = periodCell{i}(:,end);
    if i==1
        modeAmpl = amplCell{i}(:,end)/amplCell{i}(end,end);
    else
        modeAmpl = amplCell{i}(:,end)/amplCell{i}(1,end);
    end
    if size(periodCell{i},1)==1
        plot(testsCell{i}, avgPeriod/3600, 'Color', 'b', 'LineWidth', 1)
    else
        prctilePlot(testsCell{i}, periodCell{i}(:,1:end-1)/3600)
        hold on
        % plot(testsCell{i}, modePeriod/3600)
        if spaghetti
            for k = 1:4:size(periodCell{i},2) %#ok<UNRCH>
                plot(testsCell{i}, periodCell{i}(:,k)/3600,'LineWidth',.2,'Color',[0,0,1,0.1])
            end
        end
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
        legend('Estimated parameter distribution','','', '','', 'Experiments')
    end
    ylim([20 26.5])
    prettyGraph
    subplot(2,5,5+i);
    if size(periodCell{i},1)==1
        plot(testsCell{i}, avgAmpl/normAmpl, 'Color', 'b', 'LineWidth', 1)
    else
        prctilePlot(testsCell{i}, amplCell{i}(:,1:end-1)/normAmpl)
        hold on
        % plot(testsCell{i}, modeAmpl)
        if spaghetti
            for k = 1:4:size(amplCell{i},2) %#ok<UNRCH>
                plot(testsCell{i}, amplCell{i}(:,k)/normAmpl,'LineWidth',.2,'Color',[0,0,1,0.1])
            end
        end
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
    ylim([0 3.2])
    % set(gca,'YScale','log')
end

%% period and ampl summary with spaghetti plots
figure
for i = 1:length(testsCell)
    % plot periods
    subplot(2,5,i)
    hold on
    for k = 1:size(periodCell{i},2)
        plot(testsCell{i}, periodCell{i}(:,k)/3600,'Color',[0,0,1,0.1],'LineWidth',.2)
    end
    if i==1
        set(gca,'XScale','log')
        xlim([min(testsCell{i}), max(testsCell{i})])
        refIdx = length(testsCell{i});
    else
        xlim([-.05*max(testsCell{i}), max(testsCell{i})])
        refIdx = 1;
    end
    errorbar(expCasesCell{i}, expPeriodCell{i}(:,1),...
        expPeriodCell{i}(:,2), expPeriodCell{i}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r','LineWidth',1.5)
    xlabel(xlabelCell{i})
    ylabel('Period (hr)')
    ylim([20 27])
    prettyGraph
    % plot ampl
    subplot(2,5,i+5)
    hold on
    for k = 1:size(amplCell{i},2)
        plot(testsCell{i}, amplCell{i}(:,k)/amplCell{i}(refIdx,k),...
            'Color',[0,0,1,0.1],'LineWidth',.2)
    end
    if i==1
        set(gca,'XScale','log')
        xlim([min(testsCell{i}), max(testsCell{i})])
    else
        xlim([-.05*max(testsCell{i}), max(testsCell{i})])
    end
    errorbar(expCasesCell{i}, expAmplCell{i}(:,1),...
        expAmplCell{i}(:,2), expAmplCell{i}(:,2),...
        'LineStyle','none','LineWidth',1,'Marker','s','Color','r','LineWidth',1.5)
    xlabel(xlabelCell{i})
    ylabel('Relative amplitude')
    ylim([0 4])
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
                 'Control', '10μM Y27632', '20μM Y27632',...
                 'Control', '1μM CytD', '2μM CytD', '5μM CytD',...
                 'Control', '1μM LatB', '2μM LatB',...
                 'Control', '0.1μM Jas', '0.2μM Jas', '0.5μM Jas'};
groups = {1:3, 4:6, 7:10, 11:13, 14:17};
for k = 1:length(groups)
    figure
    hold on
    prettyGraph
    for i = groups{k}
        ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/pSol(18)));
        kra_cur = 1 / (1 + (LatBTests(i)/pSol(26))) + (1 + pSol(27))*JasTests(i) / pSol(28);
        inhibVec = [kra_cur, ROCKLevel, 1, 0, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD]
        inhibVec(6:9) = [1,1,3000,1];
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
gof = goodnessOfFit(periodTest', periodVec', 'NMSE');

function logLikelihood = mechanoCircadian_logLikelihoodDyn(p, y)
    % compute log likelihood for given set of parameters (stored in p
    % matrix - each row represents a separate set of parameters)
    % note that parfor is used to speed up evaluation
    % dynamics data is stored in y.dataVec
    % order of tests: diff stiffness, ROCK inhibitor, CytoD, Latrunculin B, Jas
    data = y;
    p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
        0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
        100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
    fixedParam = [5, 14, 17, 22, 25, 44:45];
    varyLogic = true(length(p0),1);
    varyLogic(fixedParam) = false;
    logLikelihood = zeros(size(p,1),1);
    parfor k = 1:size(p,1)
        pCur = p0;
        pCur(varyLogic) = p(k,:).*pCur(varyLogic)';
        logLikelihood(k) = computeLogLikelihood(pCur, data);
    end
end

function LL = computeLogLikelihood(pCur, data, varargin)
    dataVec = data.dataVec;
    periodVec = data.periodVec;
    periodStdVec = data.periodStdVec;
    if length(varargin) == 1
        plotLogic = varargin{1};
    else
        plotLogic = false;
    end
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
    expMeanDyn = dataVec(1:length(JasTests), :); %#ok<PFBNS>
    expStdDyn = dataVec(length(JasTests)+1:end, :);
    amplTest = zeros(size(JasTests));

    modelCases = [1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17]; % only non-repeats in conditions
    modelStored = cell(size(JasTests));
    minVals = zeros(size(JasTests));
    periodVals = zeros(size(JasTests));
    for i = modelCases
        ROCKLevel =  1 / (1 + (ROCKInhibTests(i)/pCur(18)));
        kra_cur = 1 / (1 + (LatBTests(i)/pCur(26))) + (1 + pCur(27))*JasTests(i) / pCur(28);
        inhibVec = [kra_cur, ROCKLevel, 1, 0, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD]
        inhibVec(6:9) = [1,1,3000,1];
        [periodCur, amplCur, t, y] = conditionToOutputs(pCur, stiffnessTests(i), inhibVec);
        % if i == 3
        %     testIdx = find(t > 60*60*24);
        %     relControlAmplB = range(y(testIdx,1))/min(y(testIdx,1));%amplCur(1)/min(y(:,1));
        %     relControlAmplP = range(y(testIdx,2))/min(y(testIdx,2));%amplCur(2)/min(y(:,2));
        % end
        amplTest(i) = amplCur(4);
        periodVals(i) = periodCur(4);
        oscDynamics = y(:,4);
        modelStored{i} = [t, oscDynamics];
        minVals(i) = min(oscDynamics);
    end
    amplTest([4,7,11,14]) = amplTest(1); % all control cases
    minVals([4,7,11,14]) = minVals(1);
    periodVals([4,7,11,14]) = periodVals(1);
    amplTestRaw = amplTest;
    % amplTest = amplTest/amplTest(1); % normalize to control amplitude
    modelStored{4} = modelStored{1};
    modelStored{7} = modelStored{1};
    modelStored{11} = modelStored{1};
    modelStored{14} = modelStored{1};

    LLvec = zeros(size(JasTests));
    % nVec = [4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
    % figure
    if plotLogic
        figure
    end
    sumAll = 0;
    lenAll = 0;
    for i = 1:length(JasTests)
        sumAll = sumAll + sum(modelStored{i}(:,2)); 
        lenAll = lenAll + length(modelStored{i}(:,2));
    end
    meanAll = sumAll/lenAll;
    for i = 1:length(JasTests)
        yNormalized = (modelStored{i}(:,2) - mean(modelStored{i}(:,2)))/amplTestRaw(1);
        % yNormalized = (modelStored{i}(:,2) - min(minVals))/amplTestRaw(1);
        % yNormalized = (modelStored{i}(:,2) - min(modelStored{i}(:,2)))/amplTestRaw(i);
        tTest = modelStored{i}(:,1);
        curMeanDyn = expMeanDyn(i,1:length(tTest))';
        curStdDyn = expStdDyn(i,1:length(tTest))';
        % curMeanDyn = curMeanDyn - min(curMeanDyn-curStdDyn);
        useLogic = curStdDyn > 0; % in certain cases, this may be zero
        LLvec(i) = -sum(((yNormalized(useLogic) - curMeanDyn(useLogic))./curStdDyn(useLogic)).^2);
        if plotLogic
            subplot(3,6,i)
            plot(tTest,curMeanDyn)
            hold on
            plot(tTest,yNormalized)
            title(sprintf("%d",i))
        end
    end
    LL = sum(LLvec);
    weightPeriod = 100;
    LL = LL - weightPeriod*sum(((periodVals/3600-periodVec)./periodStdVec).^2);
    % if relControlAmplB > 0.2 || relControlAmplP < 0.5
    %     LL = -1e7;
    % else
    % LL = LL - 1e6*(relControlAmplB-0.1)^2 - ...
    %     1e6*(relControlAmplP-1.0)^2;
    % end
end

function [periodTests, amplTests] = plotConditions(paramMat, testsMat, plotDynMat, expDataCell, legendEntries)
% testsMat has numTests columns and 6 rows: stiffness, ROCK inhib, CytD, LatB, Jas
% plotDynMat has similar structure, but contains entries to plot (number of
% curves given by number of columns in the matrix)
    
    numTests = size(testsMat,2);
    numSamples = size(paramMat,1);
    periodTests = zeros(numTests,numSamples);
    amplTests = zeros(numTests,numSamples);
    for k = 1:numSamples
        pCur= paramMat(k,:);
        for i = 1:numTests
            ROCKLevel =  1 / (1 + (testsMat(2,i)/pCur(18)));
            kra_cur = 1 / (1 + (testsMat(4,i)/pCur(26))) +...
                (1 + pCur(27))*testsMat(5,i)/pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 0, testsMat(3,i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
            inhibVec(6:9) = [1,1,3000,1];
            [curPeriod, curAmpl,~,~,~,curOscDecay] = conditionToOutputs(pCur, testsMat(1,i), inhibVec, 3600*24*7);
            if curOscDecay(4) < -.02
                if i > 1
                    periodTests(i,k) = periodTests(i-1,k);
                    amplTests(i,k) = amplTests(i-1,k);
                else
                    periodTests(i,k) = nan;
                    amplTests(i,k) = nan;
                end
            else
                periodTests(i,k) = curPeriod(4);
                amplTests(i,k) = curAmpl(4);
            end
        end
    end

    figure
    sp3 = subplot(1,3,3);
    xlabel('Time (days)')
    ylabel('PER:Luc luminescence (au)')
    prettyGraph
    hold on
    oscDynamics = cell(size(plotDynMat,2),1);
    % figure
    for k = 1:size(paramMat,1)
        pCur= paramMat(k,:);
        for i = 1:size(plotDynMat,2)
            ROCKLevel =  1 / (1 + (plotDynMat(2,i)/pCur(18)));
            kra_cur = 1 / (1 + (plotDynMat(4,i)/pCur(26))) +...
                (1+pCur(27))*plotDynMat(5,i) / pCur(28);
            inhibVec = [kra_cur, ROCKLevel, 1, 0, plotDynMat(3,i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
            inhibVec(6:9) = [1,1,3000,1];
            [~, ~, t, y] = conditionToOutputs(pCur, plotDynMat(1,i), inhibVec, 3600*10*24);
            storeLogic = t <= 24*5*3600;
            if max(t) < 24*5*3600
                oscDynamics{i}(k,:) = nan;
            else
                oscDynamics{i}(k,:) = y(storeLogic,4);
            end
            % plot(t,y(:,3))
            % hold on
            % drawnow
        end
    end
    tStored = t(storeLogic);
    meanDyn = zeros(size(plotDynMat,2),length(tStored));
    stdDyn = zeros(size(plotDynMat,2),length(tStored));
    popLogic = size(paramMat,1) > 1;
    % colororder = linspecer(size(plotDynMat,2));
    for i = 1:size(plotDynMat,2)
        if ~popLogic
            plot(tStored/(24*3600), oscDynamics{i},'LineWidth',1)
        else
            meanDyn(i,:) = mean(oscDynamics{i}(1:end-1,:),1);
            stdDyn(i,:) = std(oscDynamics{i}(1:end-1,:),1)/sqrt(size(oscDynamics{i},1));
            % errorbar(tStored/(24*3600), meanDyn(i,:), stdDyn(i,:), stdDyn(i,:))
            prctilePlot(tStored'/(24*3600), oscDynamics{i}(1:end-1,:)')
            % for k = 1:size(oscDynamics{i},1)-1
            %     plot(tStored/(24*3600), oscDynamics{i}(k,:),'Color',[colororder(i,:),0.1])
            % end
        end
    end
    set(gca,'YScale','log')
    xlim([0 5])
    legend(legendEntries)

    sp1 = subplot(1,3,1);
    activeTestsLogic = false(size(testsMat,1),1);
    for i = 1:size(testsMat,1)
        activeTestsLogic(i) = ~all(testsMat(i,:)==testsMat(i,1)); 
    end
    testStrings = {'Substrate stiffness (kPa)', 'Y27632 concentration (μM)'...
        'CytD concentration (μM)', 'LatB concentration (μM)', 'Jasplakinolide concentration (μM)'};
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
        plot(testsVec, periodTests/3600, 'Color', 'b', 'LineWidth', 1)
    else
        prctilePlot(testsVec, periodTests(:,1:end-1)/3600)
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
        plot(testsVec, amplTests/amplTests(controlIdx), 'Color', 'b', 'LineWidth', 1)
    else
        prctilePlot(testsVec, amplTests(:,1:end-1)/median(amplTests(controlIdx,1:end-1)))
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