%% Test mechano-Circadian model for range of stiffnesses
stiffnessVals = logspace(-1,3,20);
inhibMag = 0; % optionally set inhibitor case (currently, this is set to change Jas concentration)
[stiffnessMesh,inhibMesh] = meshgrid(stiffnessVals,inhibMag);
stiffnessVals = stiffnessMesh(:);
inhibMag = inhibMesh(:);
maxTime = 3600*480;
period = zeros(size(stiffnessVals));
amplitude = zeros(size(stiffnessVals));
oscDecayRate = zeros(size(stiffnessVals));
YAPTAZEq = zeros(size(stiffnessVals));
MRTFEq = zeros(size(stiffnessVals));
FcytoEq = zeros(size(stiffnessVals));
GactinEq = zeros(size(stiffnessVals));
figure
hold on
colorSeries = [0,0,.8; .796,0,.8; 0,.69,.314; 1,0,0]; % color scheme from plots
noiseLevel = [0,0]; %optional - test noise for SDDE solver

for i = 1:length(stiffnessVals)
    latBConc = 0;
    jaspConc = inhibMag(i);
    actinInhib =  1 / (1 + (latBConc/pSol(26))) + (1 + pSol(27))*jaspConc / pSol(28);
    ROCKInhib = 1;
    MRTFInhib = 1;
    YAPOverexpress = 0;
    cytoDConc = 0;
    LATSFactor = 1;
    inhibVec = [actinInhib, ROCKInhib, MRTFInhib, YAPOverexpress, cytoDConc, LATSFactor]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)]
    inhibVec(7:9) = [1,3000,1];
    [periodCur, amplitudeCur, ~, ~, rawOutput,oscDecayRateCur] = conditionToOutputs(pSol,stiffnessVals(i),inhibVec,maxTime,0,noiseLevel);
    period(i) = periodCur(3);
    amplitude(i) = amplitudeCur(3);
    oscDecayRate(i) = oscDecayRateCur(3);
    T = rawOutput{1};
    Y = rawOutput{2};
    ySS = rawOutput{3};
    oscVarIdx = 2;
    if mod(i-1,1)==0 % plot certain cases
        refConc = 13.5;
        oscDynamics = Y(:,oscVarIdx);
        [~,locs] = findpeaks(oscDynamics);
        tShift = T(locs(2))/(24*3600);
        TInterp = 0:15*60:maxTime;
        oscDynamics = interp1(T, oscDynamics, TInterp');
        plot(TInterp/(24*3600) - tShift,refConc*oscDynamics,'LineWidth',1,...
            'LineStyle','-','Marker','none')%,'Color',colorSeries(i,:))
        xlim([0 5])
        ylabel('Nuclear concentration (nM)')
        xlabel('Time (days)')
        prettyGraph
    end
    eqTimeIdx = find(T>72*3600);
    YAPTAZEq(i) = ySS(15)/(ySS(17)+ySS(18));
    MRTFEq(i) = ySS(25)/ySS(26);
    FcytoEq(i) = ySS(5);
    GactinEq(i) = ySS(9);
end

%% plot detailed control case (Fig 1C)
refConc = 13.5;
maxTime = 240*3600;
inhibVec = [1, 1, 1, 1, 0, 1]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), cytD, LATS] 
[t, y, ySS] = MechanoCircadianModel([0 maxTime], [1e7,inf], pSol, inhibVec, 0);%kraMult(i));
[~,locs] = findpeaks(y(:,2));
tShift = t(locs(3))/(3600);
figure
subplot(3,1,1)
plot(t/3600 - tShift,refConc*y(:,1), 'LineWidth', 1.5)
ylim([0 0.5*refConc])
xlabel('Time (hr)')
ylabel({'Nuclear BMAL1', '(nM)'})
hold on
yyaxis right
plot(t/3600 - tShift,refConc*y(:,2), 'LineWidth', 1.5)
ylim([0 1.5*refConc])
ylabel({'Nuclear PER/CRY', '(nM)'})
xlim([0 48])
prettyGraph
subplot(3,1,2)
KeB = pSol(3)./(1 + (y(:,1)./pSol(4)).^pSol(2));
YConc = ySS(15);
MConc = ySS(25);
KeB2 = (pSol(12)*YConc^2/(pSol(13)^2+YConc^2) + pSol(23)*MConc^2/(pSol(24)^2+MConc^2));
KeP = pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
KeP2 = (pSol(15)*MConc^2/(pSol(16)^2+MConc^2) + pSol(20)*YConc^2/(pSol(21)^2+YConc^2));
plot((t+pSol(1))/3600 - tShift, refConc*(KeB+KeB2)*3600, 'LineWidth', 1.5)
xlabel('Time (hr)')
ylabel({'BMAL1 expression','rate (nM/hr)'})
ylim([0 0.5*refConc])
hold on
yyaxis right
plot((t+pSol(7))/3600 - tShift, refConc*(KeP+KeP2)*3600, 'LineWidth', 1.5)
ylim([0 2.5*refConc])
ylabel({'PER/CRY expression','rate (nM/hr)'})
xlim([0 48])
prettyGraph
subplot(3,1,3)
KdBP = refConc*3600*pSol(5)*y(:,1).*y(:,2);
plot(t/3600 - tShift, KdBP, 'r', 'LineWidth', 1.5)
ylim([0 0.35*refConc])
xlim([0 48])
xlabel('Time (hr)')
ylabel({'B-P degradation', 'rate (nM/hr)'})
prettyGraph

figure
plot(t/(3600*24) - tShift/24,refConc*y(:,1), 'LineWidth', 1.5)
ylim([0 0.45*refConc])
xlabel('Time (hr)')
ylabel({'Nuclear BMAL1 (nM)'})
yyaxis right
plot(t/(3600*24) - tShift/24,refConc*y(:,2), 'LineWidth', 1.5)
ylim([0 1.75*refConc])
ylabel({'Nuclear PER/CRY (nM)'})
prettyGraph
xlim([0 5])

%% generate populations for Figs 4-5
fig4 = false;
fig5 = true;
if fig5
    stiffnessVals = [30, 0.3, 30, 30, 30, 30, 30, 30, 30];
    cytDConc = [0, 0, 1, 0, 0, 0, 0, 0, 0];
    latAConc = [0, 0, 0, 0.2, 0, 0, 0, 0, 0];
    LATSFactor = [1, 1, 1, 1, 7.5, 1, 1, 1, 1]; % 1 is low density, 7.5 for high density
    blebbiConc = [0, 0, 0, 0, 0, 10, 0, 0, 0];
    jaspConc = [0, 0, 0, 0, 0, 0, 1, 0, 0];
    contactArea = [3000, 1000, 5000, 600, 1200, 4000, 3000, 1600, 900];
elseif fig4 %#ok<UNRCH>
    stiffnessVals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1e7];
    cytDConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    latAConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    LATSFactor = [1, 1, 1, 1, 1, 1, 1, 1, 1]; % 1 is low density, 7.5 for high density
    blebbiConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    jaspConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    contactArea = 3000*ones(1,9);
end
popVar = .04; % started by testing 0.04 variance
numCells = 500;
samples = myBayesianAnalysis.Results.Sample(501:end,:,:);
numParam = size(samples,2);
popParam = ones(numCells, 34);
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 1/600;...
    100; 1; 10; 2; 2];
fixedParam = [2, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
flatSamples = zeros(numParam, size(samples,1)*size(samples,3));
for i = 1:size(samples,1)
    flatSamples(:,size(samples,3)*(i-1)+1:size(samples,3)*i) = samples(i,:,:);
end
rVals = randi(size(flatSamples,2), numCells, 1);
for i = 1:numCells
    popParam(i,varyLogic) = flatSamples(:,rVals(i));
end
popParam = popParam .* p0';
randPop = sqrt(popVar)*randn([numCells, 112]);

maxTime = 3600*1000;
Fs = 1/(15*60); % match experiment case of measuring every 15 minutes
tInterp = 0:1/(3600*Fs):960;
N = length(tInterp); 
freq = 0:Fs/N:Fs/2;
noiseLevel = [0,0,0]; %optional
snr = 5;
oscVarIdx = 1; % look at BMAL1 here

popSeq = cell(size(stiffnessVals));
meanPeriod = zeros(size(stiffnessVals));
stdPeriod = zeros(size(stiffnessVals));
percentOsc = zeros(size(stiffnessVals));
figure
hold on
for k = 1:length(stiffnessVals)
    oscStored = zeros(numCells,length(tInterp));
    allPower = zeros(numCells, floor(N/2+1));
    YAPTAZStored = zeros(numCells,1);
    MRTFStored = zeros(numCells,1);
    period = zeros(numCells,1);
    ampl = zeros(numCells,1);
    oscDecayRate = zeros(numCells,1);
    oscLogic = true(numCells,1);
    for i = 1:numCells
        pCur = popParam(i,:);
        inhibVec = ones(1,9);
        inhibVec(1) =  1 / (1 + (latAConc(k)/popParam(i,26))) + (1 + popParam(i,27))*jaspConc(k) / popParam(i,28);
        inhibVec(4) = 0*exp(randPop(i,112));%0.1; % fold overexpression of 5SA-YAP
        inhibVec(5) = cytDConc(k); % CytD
        inhibVec(6) = LATSFactor(k)*exp(randPop(i,109)); %high density LATS factor=7.5, low density=1
        blebbSens = 1.0*exp(randPop(i,110));
        inhibVec(7) = 1/(1 + blebbiConc(k)/blebbSens);
        inhibVec(8) = contactArea(k)*exp(randPop(i,111));
        inhibVec(9) = 0; %0% LaminA phosphorylation
        [periodCur,amplCur,tCur,yCur,rawOutput,oscDecayCur] =...
            conditionToOutputs(pCur,stiffnessVals(k),inhibVec,maxTime,randPop(i,1:105),noiseLevel);
        period(i) = periodCur(oscVarIdx);
        ampl(i) = amplCur(oscVarIdx);
        oscDecayRate(i) = oscDecayCur(oscVarIdx);
        if oscDecayCur(oscVarIdx) < -8e-4
            oscLogic(i) = false;
        end
        T = rawOutput{1};
        Y = rawOutput{2};
        ySS = rawOutput{3};
        oscCur = interp1(T/3600, Y(:,oscVarIdx), tInterp);
        oscCur = awgn(oscCur, snr); % add white Gaussian noise
        oscCur(oscCur<0) = 0;
        if ~mod(i,10)
            fprintf('Cell %d of %d, test %d\n', i, numCells, k)
        end
        oscStored(i,:) = oscCur;
        YAPTAZStored(i) = ySS(15)/(ySS(17)+ySS(18));
        MRTFStored(i) = ySS(25)/ySS(26);
        yPowerCalc = oscCur';
        yPowerCalc = filtfilt([0.2,0.2,0.2,0.2,0.2], 1, yPowerCalc); % same filter settings as Abenza et al
        yPowerCalc = (yPowerCalc - mean(yPowerCalc)) / std(yPowerCalc);
        Y_fft = fft(yPowerCalc);
        Y_fft = Y_fft(1:floor(N/2+1));
        curPower = (1/(Fs*N)) * abs(Y_fft).^2;
        curPower(2:end-1) = 2*curPower(2:end-1);
        allPower(i,:) = curPower;
    end
    meanPeriod(k) = mean(period(period<maxTime)/3600);
    stdPeriod(k) = std(period(period<maxTime)/3600);
    percentOsc(k) = sum(oscLogic) / length(oscLogic);
    % Compute the "Circadian power fraction"
    avgPower = mean(allPower,1);
    freqSelLogic = freq > 0.7/(24*3600) & freq < 1.3/(24*3600);
    powerSel = avgPower(freqSelLogic);
    freqSel = freq(freqSelLogic);
    [~,maxIdx] = max(powerSel);
    startFit = max([1,maxIdx-2]);
    endFit = min([length(freqSel),maxIdx+2]);
    fitAtMax = polyfit(freqSel(startFit:endFit), powerSel(startFit:endFit), 2);
    rootFreq = roots([2*fitAtMax(1), fitAtMax(2)]);
    if rootFreq > min(freqSel) && rootFreq < max(freqSel)
        peakFreq = rootFreq;
    else
        peakFreq = freq(maxIdx);
    end
    powerFraction = zeros(size(allPower,1),1);
    for i = 1:size(allPower,1)
        totPower = trapz(freq, allPower(i,:));
        freqInterp = peakFreq - 0.2/(3600*24):0.001/(3600*24):peakFreq + 0.2/(3600*24);
        powerInterp = interp1(freq, allPower(i,:), freqInterp);
        powerFraction(i) = trapz(freqInterp, powerInterp) / totPower;
    end
    popSeq{k} = {oscStored, period, ampl, YAPTAZStored, MRTFStored, refPeriod, powerFraction, oscDecayRate, avgPower};
    ksdensity(powerFraction,'Support',[0 1],'BoundaryCorrection','reflection')
    drawnow
end

%% plot single cell distributions for Fig4 (stiffness tests)
figure
% subplot(2,1,1)
violinplot([popSeq{9}{2},popSeq{7}{2},popSeq{5}{2},popSeq{3}{2},popSeq{1}{2}]/3600)
xticklabels({'0.1','1','10','100','1e7'})
ylim([23 28.5])
xlabel('Substrate stiffness (kPa)')
ylabel('Oscillation period (hr)')
prettyGraph

figure
% subplot(2,1,2)
violinplot([popSeq{9}{7},popSeq{7}{7},popSeq{5}{7},popSeq{3}{7},popSeq{1}{7}])
xticklabels({'0.1','1','10','100','1e7'})
ylim([0 1])
xlabel('Substrate stiffness (kPa)')
ylabel('Circadian power fraction')
prettyGraph

%% anova on populations for stiffness (Fig 4)
numCells = length(popSeq{1}{2});
conditionsVec = [0.1*ones(1,numCells), 0.3*ones(1,numCells), 1*ones(1,numCells),...
    3*ones(1,numCells), 10*ones(1,numCells), 30*ones(1,numCells),...
    100*ones(1,numCells), 300*ones(1,numCells), 1e7*ones(1,numCells)];
periodVec = [popSeq{9}{2}', popSeq{8}{2}', popSeq{7}{2}',...
    popSeq{6}{2}', popSeq{5}{2}', popSeq{4}{2}',...
    popSeq{3}{2}', popSeq{2}{2}', popSeq{1}{2}']/3600; 
powerVec = [popSeq{9}{7}', popSeq{8}{7}', popSeq{7}{7}',...
    popSeq{6}{7}', popSeq{5}{7}', popSeq{4}{7}',...
    popSeq{3}{7}', popSeq{2}{7}', popSeq{1}{7}']; 
aov_period = anova(conditionsVec, periodVec) %#ok<*NOPTS>
m_period = multcompare(aov_period)
aov_power = anova(conditionsVec, powerVec)
m_power = multcompare(aov_power) %#ok<NASGU>


%% visualize distributions at each
medianPeriod = zeros(size(stiffnessVals));
lowerBarPeriod = zeros(size(stiffnessVals));
upperBarPeriod = zeros(size(stiffnessVals));
for i = 1:length(stiffnessVals)
    distOutput = prctile(popSeq{i}{2}, [25,50,75])/3600;
    medianPeriod(i) = distOutput(2);
    lowerBarPeriod(i) = distOutput(2) - distOutput(1);
    upperBarPeriod(i) = distOutput(3) - distOutput(2);
end
figure
errorbar(latAConc, medianPeriod, lowerBarPeriod, upperBarPeriod)

%% Mimic Abenza plots
meanPowerFraction = zeros(1,length(popSeq));
stdPowerFraction = zeros(1,length(popSeq));
meanYAPTAZ = zeros(1,length(popSeq));
stdYAPTAZ = zeros(1,length(popSeq));
meanMRTF = zeros(1,length(popSeq));
stdMRTF = zeros(1,length(popSeq));
allYAPTAZ = [];
allPowerFraction = [];
allMRTF = [];
MRTFFig = figure;
title('MRTF/PowerFraction correlation')
prettyGraph
hold on
YAPTAZFig = figure;
title('YAPTAZ/PowerFraction correlation')
prettyGraph
hold on

conditionsStrings = {'Control','Soft (0.3 kPa) substrate','1 uM CytD', '0.2 uM LatA',...
    'High cell density', '10 uM blebbistatin', '1 uM Jasplakinolide', '1600 um2 micropatterns', '900 um2 micropatterns'};

inclIdx = [1:6,8,9];
colorOrder = linspecer(length(popSeq),'Qualitative');
for i = inclIdx
    curPowerFraction = popSeq{i}{7};
    curMRTF = popSeq{i}{5};
    curYAPTAZ = popSeq{i}{4};
    if i==2
        curPowerFraction = popSeq{1}{7};
    end
    meanPowerFraction(i) = median(curPowerFraction);
    stdPowerFraction(i) = std(curPowerFraction)/sqrt(length(curPowerFraction));
    meanYAPTAZ(i) = median(curYAPTAZ);
    stdYAPTAZ(i) = std(curYAPTAZ)/sqrt(length(curYAPTAZ));
    meanMRTF(i) = median(curMRTF);
    stdMRTF(i) = std(curMRTF)/sqrt(length(curMRTF));
    allYAPTAZ = [allYAPTAZ; popSeq{i}{4}]; %#ok<AGROW>
    allMRTF = [allMRTF; popSeq{i}{5}]; %#ok<AGROW>
    allPowerFraction = [allPowerFraction; popSeq{i}{7}]; %#ok<AGROW>
    figure(YAPTAZFig)
    YAPTAZQuart = prctile(curYAPTAZ,[40,50,60]);
    MRTFQuart = prctile(curMRTF,[40,50,60]);
    powerFractionQuart = prctile(curPowerFraction,[40,50,60]);
    figure(YAPTAZFig)
    errorbar(YAPTAZQuart(2),powerFractionQuart(2),diff(powerFractionQuart(1:2)),diff(powerFractionQuart(2:3)),...
        diff(YAPTAZQuart(1:2)), diff(YAPTAZQuart(2:3)),'LineWidth',1,'Marker','o','LineStyle','none','Color',colorOrder(i,:))
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    figure(MRTFFig)
    errorbar(MRTFQuart(2),powerFractionQuart(2),diff(powerFractionQuart(1:2)),diff(powerFractionQuart(2:3)),...
        diff(MRTFQuart(1:2)), diff(MRTFQuart(2:3)),'LineWidth',1,'Marker','o','LineStyle','none','Color',colorOrder(i,:))
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    % figure(YAPTAZFig)
    % errorbar(meanYAPTAZ(i), meanPowerFraction(i), stdPowerFraction(i), stdPowerFraction(i),...
    %         stdYAPTAZ(i), stdYAPTAZ(i),'LineStyle','none','LineWidth',1,'Marker','o')
    % figure(MRTFFig)
    % errorbar(meanMRTF(i), meanPowerFraction(i), stdPowerFraction(i), stdPowerFraction(i),...
    %         stdMRTF(i), stdMRTF(i),'LineStyle','none','LineWidth',1,'Marker','o')
end
[rhoYAPTAZ, pYAPTAZ] = corr(log(meanYAPTAZ(inclIdx)'), log(meanPowerFraction(inclIdx)'))
polyYAPTAZ = polyfit(log(meanYAPTAZ(inclIdx)), log(meanPowerFraction(inclIdx)), 1);
figure(YAPTAZFig)
YAPTAZRange = [min(meanYAPTAZ(inclIdx)), max(meanYAPTAZ(inclIdx))];
plot(YAPTAZRange, exp(polyval(polyYAPTAZ,log(YAPTAZRange))))
xlabel('YAP/TAZ N/C')
ylabel('Circadian power fraction')

[rhoMRTF, pMRTF] = corr(log(meanMRTF(inclIdx)'), log(meanPowerFraction(inclIdx)'))
polyMRTF = polyfit(log(meanMRTF(inclIdx)), log(meanPowerFraction(inclIdx)), 1);
figure(MRTFFig)
MRTFRange = [min(meanMRTF(inclIdx)), max(meanMRTF(inclIdx))];
plot(MRTFRange, exp(polyval(polyMRTF,log(MRTFRange))))
legend(conditionsStrings(inclIdx))
xlabel('MRTF N/C')
ylabel('Circadian power fraction')

%% plot population dynamics in heat maps
figure
plotIdx = [1,2,3];
for i = 1:length(plotIdx)
    subplot(1,length(plotIdx),i)
    imagesc(popSeq{plotIdx(i)}{1}*13.5)
    colorbar
    ylim([0 200])
    xticks(1+4*24*[0, 1, 2, 3, 4, 5, 6, 7])
    xticklabels({'0','1','2','3','4','5','6','7'})
    % colormap('turbo')
    clim([0 60])
    xlim([0 4*5*24])
    prettyGraph
    ylabel('Cell number')
    xlabel('Time (days)')
end

%% Compare mutant results
figure
% subplot(3,1,1)
violinplot([popSeq{1}{4},popSeq{2}{4},popSeq{3}{4}])
xticklabels({'WT','YAP mutant','LMNA mutant'})
% ylim([-.5 .4])
ylabel('YAP/TAZ N/C')
prettyGraph

figure
% subplot(3,1,2)
violinplot([popSeq{1}{5},popSeq{2}{5},popSeq{3}{5}])
xticklabels({'WT','YAP mutant','LMNA mutant'})
% ylim([-1.5 .6])
ylabel('MRTF N/C')
prettyGraph

figure
% subplot(3,1,3)
violinplot([popSeq{1}{7},popSeq{2}{7},popSeq{3}{7}])
xticklabels({'WT','YAP mutant','LMNA mutant'})
ylim([0 1])
ylabel('Circadian power fraction')
prettyGraph

%% ANOVA and multiple comparisons for mutant results
YAPTAZMat = [popSeq{1}{4}, popSeq{2}{4}, popSeq{3}{4}]; 
MRTFMat = [popSeq{1}{5}, popSeq{2}{5}, popSeq{3}{5}]; 
powerMat = [popSeq{1}{7}, popSeq{2}{7}, popSeq{3}{7}];
% oscDecayMat = [popSeq{1}{8}, popSeq{2}{8}, popSeq{3}{8}]
aov_YAPTAZ = anova(YAPTAZMat)
m_YAPTAZ = multcompare(aov_YAPTAZ)
aov_MRTF = anova(MRTFMat)
m_MRTF = multcompare(aov_MRTF)
aov_power = anova(powerMat)
m_power = multcompare(aov_power)
% aov_oscDecay = anova(oscDecayMat)
% m_oscDecay = multcompare(aov_oscDecay)