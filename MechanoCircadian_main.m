%% Test mechano-Circadian model under different conditions
stiffnessVals = logspace(-1,3,20);
inhibMag = 0;%0:.02:1;%0:.05:5;
[stiffnessMesh,inhibMesh] = meshgrid(stiffnessVals,inhibMag);
stiffnessVals = stiffnessMesh(:);
inhibMag = inhibMesh(:);
jaspConc = inhibMag;

% kraMult = logspace(-1,0.5,50);
maxTime = 3600*480;
coupleParam = [0.05/3600, 1.5, 2;
               0.05/3600, 1.5, 2]; % Hill function parameters for YAP-TAZ -> xClock [magn, Kd, n]
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
% colorSeries = [0,0,0; 0,.25,.66; .57,.75,.98];
% colorSeries = colororder;
noiseLevel = [0,0];%[.5e-4, .5e-4];

for i = 1:length(stiffnessVals)
%     coupleParam(1,1) = ratioVals(i)*.5/3600;
%     coupleParam(2,1) = (1-ratioVals(i))*.5/3600;
    % paramList = [9.5582, 3.2924, 1.2876, 0.0499, 0.6165, 0.4856, 11.0546, 3.0365, 1.2880, 0.6255, 0.5021, coupleParam(1,:), coupleParam(2,:), 1];
    paramList = pSol;
    actinInhib = 1 + pSol(27)*jaspConc(i)/(pSol(28) + jaspConc(i));%1 / (1 + (inhibMag(i)/paramList(18)));
    ROCKInhib = 1;% / (1 + (inhibMag(i)/2));
    MRTFInhib = 1;%inhibMag(i);
    YAPInhib = 1;%inhibMag(i);
    cytoDConc = 0;%inhibMag(i);
    LATSFactor = 1;%inhibMag(i);%7.5;
    inhibVec = [actinInhib, ROCKInhib, MRTFInhib, YAPInhib, cytoDConc, LATSFactor]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)] 
%     [T,Y, ySS] = MechanoCircadianModel([0 maxTime], [stiffnessVals(i),inf,inhibMag(i)], paramList, inhibVec, 0);%kraMult(i));
    [period(i), amplitude(i), ~, ~, rawOutput,oscDecayRate(i)] = conditionToOutputs(pSol,stiffnessVals(i),inhibVec,maxTime,0,noiseLevel);
    T = rawOutput{1};
    Y = rawOutput{2};
    ySS = rawOutput{3};
    oscVarIdx = 2;
    if mod(i-1,1)==0
        % plot(T/(24*3600),Y(:,oscVarIdx),'LineWidth',0.5,'LineStyle','-','Marker','none')%,'Color',colorSeries(i,:))
        oscDynamics = Y(:,3);%3600*pSol(9)./(1 + (pSol(10)./Y(:,1)).^pSol(8));
        plot(T/(24*3600),oscDynamics,'LineWidth',0.5,'LineStyle','-','Marker','none')%,'Color',colorSeries(i,:))
        xlim([0 5])
%         ylim([0 .5])
        ylabel('PER/CRY abundance')
        xlabel('Time (days)')
%         yyaxis right
%         plot(T/(24*3600) - 10,Y(:,1),'LineWidth',1,'LineStyle','-','Marker','none')
%         ylabel('BMAL1 abundance')
% %         ylim([0 0.7])
%         legend('PER/CRY','BMAL1')
%         prettyGraph
    end
    eqTimeIdx = find(T>72*3600);
    YAPTAZEq(i) = ySS(15)/(ySS(17)+ySS(18));%mean(Y(eqTimeIdx,2));
    MRTFEq(i) = ySS(25)/ySS(26);%mean(Y(eqTimeIdx,3));
    FcytoEq(i) = ySS(5);
    GactinEq(i) = ySS(9);
end
% figure
% surf(log10(stiffnessMesh),inhibMesh,reshape(period,size(stiffnessMesh)))
% view([0 90])
% plot(stiffnessVals,period)
% figure
% b = bar(period/3600);
% b.FaceColor = 'flat';
% prettyGraph
% for i = 1:length(period)
%     b.CData(i,:) = colorSeries(i,:);
% end
%% CytD validation
FibroblastCytDPeriod = [24.2103,    0.2060; % Control
                        24.9957,    0.1803; % 1 uM CytD
                        24.9700,    0.0258; % 2 uM CytD
                        24.6352,    0.3991];% 5 uM CytD
FibroblastCytDAmpl = [1.0000,    0.1954; % Control
                      1.2874,    0.1494; % 1 uM CytD
                      1.6552,    0.3218; % 2 uM CytD
                      1.9425,    0.3908]; % 5 uM CytD

CytDTests = 0:.2:10;
CytDPeriod = zeros(size(CytDTests));
CytDAmpl = zeros(size(CytDTests));
for i = 1:length(CytDTests)
    inhibVec = [1, 1, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [CytDPeriod(i), CytDAmpl(i)] = conditionToOutputs(pSol, 1e7, inhibVec, 480*3600);
end
figure
% subplot(1,2,1)
plot(CytDTests, CytDPeriod/3600, 'LineWidth', 1)
hold on
CytDConc = [0, 1, 2, 5];
errorbar(CytDConc, FibroblastCytDPeriod(:,1),...
    FibroblastCytDPeriod(:,2), FibroblastCytDPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s')
xlim([0 6])
% subplot(1,2,2)
% plot(CytDTests, CytDAmpl/CytDAmpl(1), 'LineWidth', 1)
% hold on
% errorbar(CytDConc, FibroblastCytDAmpl(:,1),...
%     FibroblastCytDAmpl(:,2), FibroblastCytDAmpl(:,2),...
%     'LineStyle','none','LineWidth',1,'Marker','s')
% xlim([0 6])

%% LatB validation
FibroblastLatBPeriod = [24.9389,    0.1048;
                        25.4105,    0.1921;
                        25.9694,    0.2620];
FibroblastLatBAmpl = [1.0000,    0.0235;
                      1.6118,    0.1412;
                      2.3765,    0.0588];
LatBTests = 0:.1:4;
LatBPeriod = zeros(size(LatBTests));
LatBAmpl = zeros(size(LatBTests));
for i = 1:length(LatBTests)
    inhibVec = [1/(1 + LatBTests(i)/0.5), 1, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [LatBPeriod(i), LatBAmpl(i)] = conditionToOutputs(pSol, 1e7, inhibVec, 480*3600);
end
figure
% subplot(1,2,1)
plot(LatBTests, LatBPeriod/3600, 'LineWidth', 1)
hold on
LatBConc = [0, 1, 2];
errorbar(LatBConc, FibroblastLatBPeriod(:,1),...
    FibroblastLatBPeriod(:,2), FibroblastLatBPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s')
xlim([0 6])
% subplot(1,2lim(,2)
% plot(LatBTests, LatBAmpl/LatBAmpl(1), 'LineWidth', 1)
% hold on
% errorbar(LatBConc, FibroblastLatBAmpl(:,1),...
%     FibroblastLatBAmpl(:,2), FibroblastLatBAmpl(:,2),...
%     'LineStyle','none','LineWidth',1,'Marker','s')
% xlim([0 6])

%% generate population
% lower YAPTAZ in nucleus: 1.25 uM
% higher YAPTAZ in nucleus: 2 uM
stiffnessVals = [1e7,1e7,1e7];%[30, 0.3, 30, 30, 30, 30, 30, 30, 30];%logspace(-1,4,20);
cytDConc = [0,0,0];%[0, 0, 1, 0, 0, 0, 0, 0, 0];
latAConc = [0,1,2];%[0, 0, 0, 0.2, 0, 0, 0, 0, 0];
LATSFactor = [1, 1, 1, 1, 7.5, 1, 1, 1, 1]; % 1 is low density, 7.5 for high density
blebbiConc = [0, 0, 0, 0, 0, 10, 0, 0, 0];
jaspConc = [0, 0, 0, 0, 0, 0, 1, 0, 0];
contactArea = [2000,2000,2000,2000];%[2000, 500, 2000, 500, 1200, 2000, 1000, 1600, 900];

popSeq = cell(size(stiffnessVals));
meanPeriod = zeros(size(stiffnessVals));
stdPeriod = zeros(size(stiffnessVals));
percentOsc = zeros(size(stiffnessVals));
% figure
% hold on
for k = 1:length(stiffnessVals)
    fixedParam = [2, 8, 14, 17, 22]; % don't allow exponents to vary
    popVar = 0.01; % started by testing 0.04 variance
    varVec = sqrt(popVar)*ones(23,1);
    varVec(fixedParam) = 0;
    numCells = 1000;
    maxTime = 3600*1000;
    tInterp = 0:0.25:960;
    Fs = 1/(15*60); % experimentally, measure every 15 minutes
    N = length(tInterp); 
    freq = 0:Fs/N:Fs/2;
    oscStored = [];%zeros(numCells,length(tInterp));
    allPower = zeros(numCells, floor(N/2+1));
    YAPTAZStored = zeros(numCells,1);
    MRTFStored = zeros(numCells,1);
    oscVarIdx = 1; % look at BMAL1 here
    period = zeros(numCells,1);
    ampl = zeros(numCells,1);
    oscLogic = false(numCells,1);
    noiseLevel = [0,0];%[1e-4, 1e-4];
    for i = 1:numCells
        varVecCur = varVec.*randn(size(varVec));
        pCur = pSol' .* exp(varVecCur(1:22));
        %inhib = [actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), cytD, LATSFactor, blebbiFactor, contactArea]
        inhibVec = ones(1,8); 
        inhibVec(1) = 1/(1 + latAConc(k)/.5) + jaspConc(k) / (0.1 + jaspConc(k));
        inhibVec(4) = 1;%0.1; %10% phosphorylation
        inhibVec(5) = cytDConc(k); % CytD
        inhibVec(6) = LATSFactor(k); %high density LATS factor=7.5, low density=1
        inhibVec(7) = 1/(1 + blebbiConc(k)/1.0);
        inhibVec(8) = contactArea(k)*exp(varVecCur(23));
        [period(i),ampl(i),~,~,rawOutput,oscLogic(i)] = conditionToOutputs(pCur,stiffnessVals(k),inhibVec,maxTime,popVar,noiseLevel);
        T = rawOutput{1};
        Y = rawOutput{2};
        ySS = rawOutput{3};
        %     [T,Y, ySS] = MechanoCircadianModel([0 maxTime], [stiffness,inf,0], pCur, inhibVec, popVar);
%         plot(T/3600 - 48, Y(:,oscVarIdx))
        oscCur = interp1(T/3600, Y(:,oscVarIdx), tInterp);
        YAPTAZStored(i) = ySS(15)/(ySS(17)+ySS(18));
        MRTFStored(i) = ySS(25)/ySS(26);
        yPowerCalc = oscCur';%oscStored(i,:)';
        yPowerCalc = filtfilt([0.2,0.2,0.2,0.2,0.2], 1, yPowerCalc);
        yPowerCalc = (yPowerCalc - mean(yPowerCalc)) / std(yPowerCalc);
        Y_fft = fft(yPowerCalc);
        Y_fft = Y_fft(1:floor(N/2+1));
        curPower = (1/(Fs*N)) * abs(Y_fft).^2;
        curPower(2:end-1) = 2*curPower(2:end-1);
        allPower(i,:) = curPower;
    end
    % compute reference period for zero variance
    refPeriod = conditionToOutputs(pSol,stiffnessVals(k),inhibVec,maxTime);
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
        peakFreq = 1/(refPeriod);%freq(maxIdx);
    end
    powerFraction = zeros(size(allPower,1),1);
    for i = 1:size(allPower,1)
        totPower = trapz(freq, allPower(i,:));
        freqInterp = peakFreq - 0.2/(3600*24):0.001/(3600*24):peakFreq + 0.2/(3600*24);
        powerInterp = interp1(freq, allPower(i,:), freqInterp);
        powerFraction(i) = trapz(freqInterp, powerInterp) / totPower;
    end
    popSeq{k} = {oscStored, period, ampl, YAPTAZStored, MRTFStored, refPeriod, powerFraction};
end

%% Random sampling from population
numSamples = 4;
samplePeriod = cell(size(stiffnessVals));
sampleMeans = zeros(size(stiffnessVals));
sampleStds = zeros(size(stiffnessVals));
for k = 1:length(stiffnessVals)
    for j = 1:numSamples
        curPop = randsample(popSeq{k}{2}, 50);
        evalLogic = curPop < maxTime;
        samplePeriod{k}(j) = mean(curPop(evalLogic));
    end
%     samplePeriod{k} = bootstrp(numBootstraps, @mean, popSeq{k}{2}(evalLogic));
    sampleMeans(k) = mean(samplePeriod{k})/3600;
    sampleStds(k) = std(samplePeriod{k})/3600;
end
figure
errorbar(latAConc, sampleMeans, sampleStds, sampleStds)
% set(gca,'xscale','log')

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
meanPowerFraction = zeros(size(meanPeriod));
stdPowerFraction = zeros(size(meanPeriod));
meanYAPTAZ = zeros(size(meanPeriod));
stdYAPTAZ = zeros(size(meanPeriod));
meanMRTF = zeros(size(meanPeriod));
stdMRTF = zeros(size(meanPeriod));
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
for i = 1:length(popSeq)
    meanPowerFraction(i) = mean(popSeq{i}{7});
    stdPowerFraction(i) = std(popSeq{i}{7})/sqrt(length(popSeq{i}{7}));
    meanYAPTAZ(i) = mean(popSeq{i}{4});
    stdYAPTAZ(i) = std(popSeq{i}{4})/sqrt(length(popSeq{i}{4}));
    meanMRTF(i) = mean(popSeq{i}{5});
    stdMRTF(i) = std(popSeq{i}{5})/sqrt(length(popSeq{i}{5}));
    allYAPTAZ = [allYAPTAZ; popSeq{i}{4}];
    allMRTF = [allMRTF; popSeq{i}{5}];
    allPowerFraction = [allPowerFraction; popSeq{i}{7}];
    figure(YAPTAZFig)
    YAPTAZQuart = prctile(popSeq{i}{4},[40,50,60]);
    MRTFQuart = prctile(popSeq{i}{5},[40,50,60]);
    powerFractionQuart = prctile(popSeq{i}{7},[40,50,60]);
%     figure(YAPTAZFig)
%     errorbar(YAPTAZQuart(2),powerFractionQuart(2),diff(powerFractionQuart(1:2)),diff(powerFractionQuart(2:3)),...
%         diff(YAPTAZQuart(1:2)), diff(YAPTAZQuart(2:3)),'LineWidth',1,'Marker','o')
%     figure(MRTFFig)
%     errorbar(MRTFQuart(2),powerFractionQuart(2),diff(powerFractionQuart(1:2)),diff(powerFractionQuart(2:3)),...
%         diff(MRTFQuart(1:2)), diff(MRTFQuart(2:3)),'LineWidth',1,'Marker','o')
    figure(YAPTAZFig)
    errorbar(meanYAPTAZ(i), meanPowerFraction(i), stdPowerFraction(i), stdPowerFraction(i),...
            stdYAPTAZ(i), stdYAPTAZ(i),'LineStyle','none','LineWidth',1,'Marker','o')
    figure(MRTFFig)
    errorbar(meanMRTF(i), meanPowerFraction(i), stdPowerFraction(i), stdPowerFraction(i),...
            stdMRTF(i), stdMRTF(i),'LineStyle','none','LineWidth',1,'Marker','o')
end
% figure
% errorbar(meanYAPTAZ, meanPowerFraction, stdPowerFraction, stdPowerFraction,...
%         stdYAPTAZ, stdYAPTAZ,'LineStyle','none','LineWidth',1,'Marker','s')
% figure
% errorbar(meanMRTF, meanPowerFraction, stdPowerFraction, stdPowerFraction,...
%         stdMRTF, stdMRTF,'LineStyle','none','LineWidth',1,'Marker','s')