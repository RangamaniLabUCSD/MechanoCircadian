%% Test mechano-Circadian model for range of stiffnesses
stiffnessVals = 100;%1e7;%logspace(-1,3,20);
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
% figure
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
inhibVec(7:9) = [1,3000,1];
[t, y, ySS] = MechanoCircadianModel([0 maxTime], [100,inf], pSol, inhibVec, 0);%kraMult(i));
[~,locs] = findpeaks(y(:,2));
tShift = t(locs(3))/(3600);
figure
subplot(3,1,1)
plot(t/3600 - tShift,refConc*y(:,1), 'LineWidth', 1.5)
ylim([0 0.4*refConc])
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
% plot((t+pSol(1))/3600 - tShift, refConc*(KeB+KeB2)*3600, 'LineWidth', 1.5)
plot((t)/3600 - tShift, refConc*(KeB)*3600, 'LineWidth', 1.5)
xlabel('Time (hr)')
ylabel({'BMAL1 expr.','rate (nM/hr)'})
ylim([0 0.25*refConc])
hold on
yyaxis right
% plot((t+pSol(7))/3600 - tShift, refConc*(KeP+KeP2)*3600, 'LineWidth', 1.5)
plot((t)/3600 - tShift, refConc*(KeP)*3600, 'LineWidth', 1.5)
ylim([0 2.5*refConc])
ylabel({'PER/CRY expr.','rate (nM/hr)'})
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
xlabel('Time (days)')
ylabel({'Nuclear BMAL1 (nM)'})
yyaxis right
plot(t/(3600*24) - tShift/24,refConc*y(:,2), 'LineWidth', 1.5)
ylim([0 1.75*refConc])
ylabel({'Nuclear PER/CRY (nM)'})
prettyGraph
xlim([0 5])

%% generate populations for Figs 4-5
fig4 = true;
fig5 = false;
if fig5
    stiffnessVals = [30, 0.3, 30, 30, 30, 30, 30, 30, 30]; %#ok<*UNRCH>
    cytDConc = [0, 0, 1, 0, 0, 0, 0, 0, 0];
    latAConc = [0, 0, 0, 0.2, 0, 0, 0, 0, 0];
    LATSFactor = [1, 1, 1, 1, 3, 1, 1, 1, 1]; % 1 is low density, 3 for high density
    blebbiConc = [0, 0, 0, 0, 0, 10, 0, 0, 0];
    jaspConc = [0, 0, 0, 0, 0, 0, 1, 0, 0];
    contactArea = [3000, 1000, 5000, 600, 1200, 4000, 3000, 1600, 900];
elseif fig4 
    stiffnessVals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1e7];
    cytDConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    latAConc = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    LATSFactor = [1, 1, 1, 1, 1, 1, 1, 1, 1]; % 1 is low density, 3 for high density
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
% flatSamples = zeros(numParam, size(samples,1)*size(samples,3));
% for i = 1:size(samples,1)
%     flatSamples(:,size(samples,3)*(i-1)+1:size(samples,3)*i) = samples(i,:,:);
% end
% rVals = randi(size(flatSamples,2), numCells, 1);
% for i = 1:numCells
%     popParam(i,varyLogic) = flatSamples(:,rVals(i));
% end
% popParam = popParam .* p0';
% randPop = sqrt(popVar)*randn([numCells, 112]);
popParam = popParamWT;
randPop = randPopWT;

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
    ySS_stored = zeros(numCells, 26);
    period = zeros(numCells,1);
    ampl = zeros(numCells,1);
    oscDecayRate = zeros(numCells,1);
    oscLogic = true(numCells,1);
    for i = 1:numCells
        pCur = popParam(i,:);
        inhibVec = ones(1,9);
        inhibVec(1) =  1 / (1 + (latAConc(k)/popParam(i,26))) + (1 + popParam(i,27))*jaspConc(k) / popParam(i,28);
        inhibVec(4) = 0*exp(randPop(i,112)); % fold overexpression of 5SA-YAP
        inhibVec(5) = cytDConc(k); % CytD
        inhibVec(6) = LATSFactor(k)*exp(randPop(i,109)); %high density LATS factor=7.5, low density=1
        blebbSens = 1.0*exp(randPop(i,110));
        inhibVec(7) = 1/(1 + blebbiConc(k)/blebbSens);
        inhibVec(8) = contactArea(k)*exp(randPop(i,111));
        inhibVec(9) = 0; % percent LaminA phosphorylation (0 for mutant, 1 for normal)
        inhibVec(10) = 2; % fold LaminA expression (no variation here because variation occurs in LaminA conc in MechanoSS) (1 for normal, 2 for mutant)
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
        ySS_stored(i,:) = ySS;
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
    refPeriod = [];
    popSeq{k} = {oscStored, period, ampl, YAPTAZStored, MRTFStored, refPeriod, powerFraction, oscDecayRate, avgPower, popParam, randPop, ySS_stored};
    ksdensity(powerFraction,'Support',[0 1],'BoundaryCorrection','reflection')
    drawnow
end
save('StiffnessTestSeq_LMNAMutant_sharedWTVals.mat','popSeq')

%% plot single cell distributions for Fig4 (stiffness tests)
figure
% subplot(2,1,1)
violinplot([popSeq{1}{2},popSeq{3}{2},popSeq{5}{2},popSeq{7}{2},popSeq{9}{2}]/3600)
xticklabels({'0.1','1','10','100','1e7'})
ylim([23 28.5])
xlabel('Substrate stiffness (kPa)')
ylabel('Oscillation period (hr)')
prettyGraph

figure
% subplot(2,1,2)
violinplot([popSeq{1}{7},popSeq{3}{7},popSeq{5}{7},popSeq{7}{7},popSeq{9}{7}])
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
periodVec = [popSeq{1}{2}', popSeq{2}{2}', popSeq{3}{2}',...
    popSeq{4}{2}', popSeq{5}{2}', popSeq{6}{2}',...
    popSeq{7}{2}', popSeq{8}{2}', popSeq{9}{2}']/3600; 
powerVec = [popSeq{1}{7}', popSeq{2}{7}', popSeq{3}{7}',...
    popSeq{4}{7}', popSeq{5}{7}', popSeq{6}{7}',...
    popSeq{7}{7}', popSeq{8}{7}', popSeq{9}{7}']; 
[p,t,aov_period] = anova1(periodVec, conditionsVec) %#ok<*ASGLU,*NOPTS>
m_period = multcompare(aov_period)
[p,t,aov_power] = anova1(powerVec, conditionsVec)
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

%% Mimic Abenza plots (Fig 4)
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
        curPowerFraction = (curPowerFraction + popSeq{1}{7})/2;
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

%% load color map 
% (this is tailored to a single json (inferno) currently, 
% would need to be edited for loading in other color maps)
[json_file,json_path] = uigetfile('*.json','Select json file with color map');
json_file = fullfile(json_path,json_file);
load_colors = readcell(json_file,'FileType','text');
RGBPoints = cell2mat(load_colors(16:1039));
RGBPoints = reshape(RGBPoints',[4, length(RGBPoints)/4]);
dataVals = RGBPoints(1,:);
dataInterp = dataVals(1):(dataVals(end)-dataVals(1))/1000:dataVals(end);
RGBVals = RGBPoints(2:4,:);
RInterp = interp1(dataVals,RGBVals(1,:),dataInterp);
GInterp = interp1(dataVals,RGBVals(2,:),dataInterp);
BInterp = interp1(dataVals,RGBVals(3,:),dataInterp);
RGBInterp = vertcat(RInterp,GInterp,BInterp);

%% plot population dynamics in heat maps
figure
plotIdx = [1,2];
for i = 1:length(plotIdx)
    subplot(1,length(plotIdx),i)
    imagesc(popSeq{plotIdx(i)}{1}*13.5)
    colorbar
    ylim([0 200])
    xticks(1+4*24*[0, 1, 2, 3, 4, 5, 6, 7])
    xticklabels({'0','1','2','3','4','5','6','7'})
    colormap(RGBInterp(:,150:end)')
    clim([0 60])
    xlim([0 4*3.05*24])
    prettyGraph
    ylabel('Cell number')
    xlabel('Time (days)')
end

%% average population dynamics
Fs = 1/(15*60); % match experiment case of measuring every 15 minutes
tInterp = 0:1/(3600*Fs):960;
figure
hold on
for i = 1:length(popSeq)
    meanDyn = mean(popSeq{i}{1});
    plot(tInterp,meanDyn)
    popSeq{i}{10} = meanDyn;
end

%% Compare mutant results
figure
% subplot(3,1,1)
violinplot([popSeq{1}{12}(:,15),popSeq{2}{12}(:,15),popSeq{3}{12}(:,15)])
xticklabels({'WT','YAP mutant','LMNA mutant'})
% ylim([-.5 .4])
ylabel('YAP/TAZ nuclear conc. (uM)')
prettyGraph

figure
% subplot(3,1,2)
violinplot([popSeq{1}{12}(:,25),popSeq{2}{12}(:,25),popSeq{3}{12}(:,25)])
xticklabels({'WT','YAP mutant','LMNA mutant'})
% ylim([-1.5 .6])
ylabel('MRTF nuclear conc. (uM)')
prettyGraph

figure
% subplot(3,1,3)
violinplot([popSeq{1}{7},popSeq{2}{7},popSeq{3}{7}])
xticklabels({'WT','YAP mutant','LMNA mutant'})
ylim([0 1])
ylabel('Circadian power fraction')
prettyGraph

%% ANOVA and multiple comparisons for mutant results
YAPTAZMat = [popSeq{1}{12}(:,15), popSeq{2}{12}(:,15), popSeq{3}{12}(:,15)]; 
MRTFMat = [popSeq{1}{12}(:,25), popSeq{2}{12}(:,25), popSeq{3}{12}(:,25)]; 
powerMat = [popSeq{1}{7}, popSeq{2}{7}, popSeq{3}{7}];
% oscDecayMat = [popSeq{1}{8}, popSeq{2}{8}, popSeq{3}{8}]
[p,t,aov_YAPTAZ] = anova1(YAPTAZMat)
m_YAPTAZ = multcompare(aov_YAPTAZ)
[p,t,aov_MRTF] = anova1(MRTFMat)
m_MRTF = multcompare(aov_MRTF)
[p,t,aov_power] = anova1(powerMat)
m_power = multcompare(aov_power)
% aov_oscDecay = anova(oscDecayMat)
% m_oscDecay = multcompare(aov_oscDecay)

%% movie of population-level oscillations
figure
curOscMat = popSeq{2}{1}*13.5;
for i = 1:size(curOscMat,1)
    curOscMat(i,:) = smooth(curOscMat(i,:),9);
end

xCellLocs = 0:20;
yCellLocs = 0:20;
[xCellLocs, yCellLocs] = meshgrid(xCellLocs, yCellLocs);
xCellLocs = xCellLocs(:) + 0.2*randn(size(xCellLocs(:)));
yCellLocs = yCellLocs(:) + 0.2*randn(size(yCellLocs(:)));
colorMap = parula(256);
numCells = length(yCellLocs);

saveVid = true;
if saveVid
    frameRate = 20;
    path = uigetdir('','Choose where to save video');
    myVideo = VideoWriter(fullfile(path,'YAPMutantPop_softSubstrate.mp4'),'MPEG-4');
    myVideo.FrameRate = frameRate;
    open(myVideo)
end

for i = 1:4:4*24*14+1%size(curOscMat,2)
    curOscVals = curOscMat(1:numCells,i);
    colorIdx = round(255*curOscVals/60) + 1;
    colorIdx(colorIdx > 256) = 256;
    scatter(xCellLocs, yCellLocs, 20*ones(size(xCellLocs)), colorMap(colorIdx,:), "filled")
    xlim([-1 21])
    ylim([-1 21])
    prettyGraph
    title(sprintf('t = %.2f days',(i-1)*0.25/24))
    cb = colorbar;
    cb.Label.String = 'BMAL1 (nM)'
    clim([0 60])
    xticks([])
    yticks([])

    drawnow
    if saveVid
        img = print('-RGBImage','-r300');
        writeVideo(myVideo, img)
    end
end

if saveVid
    close(myVideo)
end

%% mutant comparison for single cell with MAP parameters
stiffnessVals = [300, 300, 300];
YAPOverexpress = [0, 1, 0]
laminAPhos = [1, 1, 0];
maxTime = 3600*480;
figure
hold on
for i = 1:length(stiffnessVals)
    inhibVec = [1,1,1,YAPOverexpress(i),0,1,1,3000,laminAPhos(i)];
    [periodCur, amplitudeCur, tOut, yOut] = conditionToOutputs(pSol,stiffnessVals(i),inhibVec,maxTime,0);
    plot(tOut/3600, yOut(:,2)*13.5)
end
prettyGraph

%% mutant comparison on different stiffnesses
figure
loadWT = load('StiffnessTestSeq_WT.mat');
popSeq_WT = loadWT.popSeq;
loadLaminAMutant = load('StiffnessTestSeq_LMNAMutant_sharedWTVals.mat');
popSeq_LMNA = loadLaminAMutant.popSeq;
loadYAPMutant = load('StiffnessTestSeq_YAPMutant_sharedWTVals.mat');
popSeq_YAP = loadYAPMutant.popSeq;
stiffnessVals = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1e7];
numCells = length(popSeq_WT{1}{4});
YAPTAZ_WT = zeros(numCells,length(stiffnessVals));
YAPTAZ_LMNA = zeros(numCells,length(stiffnessVals));
YAPTAZ_YAP = zeros(numCells,length(stiffnessVals));
for i = 1:length(stiffnessVals)
    YAPTAZ_WT(:,i) = popSeq_WT{i}{4};
    YAPTAZ_LMNA(:,i) = popSeq_LMNA{i}{4};
    YAPTAZ_YAP(:,i) = popSeq_YAP{i}{4};
end
errorbar(stiffnessVals, median(YAPTAZ_WT), median(YAPTAZ_WT)-quantile(YAPTAZ_WT,0.25),...
    quantile(YAPTAZ_WT,0.75)-median(YAPTAZ_WT))
hold on
errorbar(stiffnessVals, median(YAPTAZ_LMNA), median(YAPTAZ_LMNA)-quantile(YAPTAZ_LMNA,0.25),...
    quantile(YAPTAZ_LMNA,0.75)-median(YAPTAZ_LMNA))
errorbar(stiffnessVals, median(YAPTAZ_YAP), median(YAPTAZ_YAP)-quantile(YAPTAZ_YAP,0.25),...
    quantile(YAPTAZ_YAP,0.75)-median(YAPTAZ_YAP))
% prctilePlot(stiffnessVals, YAPTAZ_WT','b')
% hold on
% prctilePlot(stiffnessVals, YAPTAZ_LMNA','r')
set(gca, 'XScale', 'log')
% subplot(3,1,1)
% violinplot([popSeq{1}{4},popSeq{2}{4},popSeq{3}{4},popSeq{4}{4},popSeq{5}{4},popSeq{6}{4}])
% xticklabels({'WT 30 kPa','WT 0.3 kPa','YAP mutant 30 kPa','YAP mutant 0.3 kPa','LMNA mutant 30 kPa','LMNA mutant 0.3 kPa'})
% ylim([-.5 .4])
ylabel('YAP/TAZ N/C')
prettyGraph

figure
MRTF_WT = zeros(length(popSeq_WT{1}{5}),length(stiffnessVals));
MRTF_LMNA = zeros(length(popSeq_WT{1}{5}),length(stiffnessVals));
MRTF_YAP = zeros(length(popSeq_WT{1}{5}),length(stiffnessVals));
for i = 1:length(stiffnessVals)
    MRTF_WT(:,i) = popSeq_WT{i}{5};
    MRTF_LMNA(:,i) = popSeq_LMNA{i}{5};
    MRTF_YAP(:,i) = popSeq_YAP{i}{5};
end
errorbar(stiffnessVals, median(MRTF_WT), median(MRTF_WT)-quantile(MRTF_WT,0.25),...
    quantile(MRTF_WT,0.75)-median(MRTF_WT))
hold on
errorbar(stiffnessVals, median(MRTF_LMNA), median(MRTF_LMNA)-quantile(MRTF_LMNA,0.25),...
    quantile(MRTF_LMNA,0.75)-median(MRTF_LMNA))
errorbar(stiffnessVals, median(MRTF_YAP), median(MRTF_YAP)-quantile(MRTF_YAP,0.25),...
    quantile(MRTF_YAP,0.75)-median(MRTF_YAP))
% prctilePlot(stiffnessVals, MRTF_WT', 'b')
% hold on
% prctilePlot(stiffnessVals, MRTF_LMNA', 'r')
set(gca, 'XScale', 'log')
% subplot(3,1,2)
% violinplot([popSeq{1}{5},popSeq{2}{5},popSeq{3}{5},popSeq{4}{5},popSeq{5}{5},popSeq{6}{5}])
% xticklabels({'WT 30 kPa','WT 0.3 kPa','YAP mutant 30 kPa','YAP mutant 0.3 kPa','LMNA mutant 30 kPa','LMNA mutant 0.3 kPa'})
% ylim([-1.5 .6])
ylabel('MRTF N/C')
prettyGraph

figure
powerFrac_WT = zeros(length(popSeq_WT{1}{7}),length(stiffnessVals));
powerFrac_LMNA = zeros(length(popSeq_WT{1}{7}),length(stiffnessVals));
powerFrac_YAP = zeros(length(popSeq_WT{1}{7}),length(stiffnessVals));
for i = 1:length(stiffnessVals)
    powerFrac_WT(:,i) = popSeq_WT{i}{7};
    powerFrac_LMNA(:,i) = popSeq_LMNA{i}{7};
    powerFrac_YAP(:,i) = popSeq_YAP{i}{7};
end
errorbar(stiffnessVals, median(powerFrac_WT), median(powerFrac_WT)-quantile(powerFrac_WT,0.25),...
    quantile(powerFrac_WT,0.75)-median(powerFrac_WT))
hold on
errorbar(stiffnessVals, median(powerFrac_LMNA), median(powerFrac_LMNA)-quantile(powerFrac_LMNA,0.25),...
    quantile(powerFrac_LMNA,0.75)-median(powerFrac_LMNA))
errorbar(stiffnessVals, median(powerFrac_YAP), median(powerFrac_YAP)-quantile(powerFrac_YAP,0.25),...
    quantile(powerFrac_YAP,0.75)-median(powerFrac_YAP))
% prctilePlot(stiffnessVals, powerFrac_WT', 'b')
% hold on
% prctilePlot(stiffnessVals, powerFrac_LMNA', 'r')
set(gca, 'XScale', 'log')
% subplot(3,1,3)
% violinplot([popSeq{1}{7},popSeq{2}{7},popSeq{3}{7},popSeq{4}{7},popSeq{5}{7},popSeq{6}{7}])
% xticklabels({'WT 30 kPa','WT 0.3 kPa','YAP mutant 30 kPa','YAP mutant 0.3 kPa','LMNA mutant 30 kPa','LMNA mutant 0.3 kPa'})
ylim([0 1])
ylabel('Circadian power fraction')
prettyGraph

%%
figure
powerMat = [powerFrac_WT(:,6), powerFrac_WT(:,2), powerFrac_YAP(:,6),...
    powerFrac_YAP(:,2), powerFrac_LMNA(:,6), powerFrac_LMNA(:,4)]; 
violinplot(powerMat(:,[1,3,5,4,6]))
xticklabels({'WT (30 kPa)','YAP mutant (30 kPa)','LMNA mutant (30 kPa)',...
    'YAP mutant (0.3 kPa)', 'LMNA mutant (3 kPa)'})
ylim([0 1])
ylabel('Circadian power fraction')
prettyGraph

%% ANOVA and multiple comparisons for mutant results
YAPTAZMat = [YAPTAZ_WT(:,6), YAPTAZ_WT(:,2), YAPTAZ_YAP(:,6),...
    YAPTAZ_YAP(:,2), YAPTAZ_LMNA(:,6), YAPTAZ_LMNA(:,4)]; 
MRTFMat = [MRTF_WT(:,6), MRTF_WT(:,2), MRTF_YAP(:,6),...
    MRTF_YAP(:,2), MRTF_LMNA(:,6), MRTF_LMNA(:,4)];  
% oscDecayMat = [popSeq{1}{8}, popSeq{2}{8}, popSeq{3}{8}]
[p,t,aov_YAPTAZ] = anova1(YAPTAZMat)
m_YAPTAZ = multcompare(aov_YAPTAZ)
[p,t,aov_MRTF] = anova1(MRTFMat)
m_MRTF = multcompare(aov_MRTF)
[p,t,aov_power] = anova1(powerMat)
m_power = multcompare(aov_power)
% aov_oscDecay = anova(oscDecayMat)
% m_oscDecay = multcompare(aov_oscDecay)