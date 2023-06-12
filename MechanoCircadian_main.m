%% YAP TAZ analysis
stiffnessVals = 1e5;%logspace(-1,3,50);
inhibMag = 0:.1:10;%0:.1:5;%0:.1:1;
[stiffnessMesh,inhibMesh] = meshgrid(stiffnessVals,inhibMag);
stiffnessVals = stiffnessMesh(:);
inhibMag = inhibMesh(:);

% kraMult = logspace(-1,0.5,50);
maxTime = 3600*2000;
coupleParam = [0.05/3600, 1.5, 2;
               0.05/3600, 1.5, 2]; % Hill function parameters for YAP-TAZ -> xClock [magn, Kd, n]
period = zeros(size(stiffnessVals));
amplitude = zeros(size(stiffnessVals));
YAPTAZEq = zeros(size(stiffnessVals));
MRTFEq = zeros(size(stiffnessVals));
FcytoEq = zeros(size(stiffnessVals));
GactinEq = zeros(size(stiffnessVals));
figure
hold on
colorSeries = [0,0,.8; .796,0,.8; 0,.69,.314; 1,0,0];
% colorSeries = [0,0,0; 0,.25,.66; .57,.75,.98];
% colorSeries = colororder;
for i = 1:length(stiffnessVals)
%     coupleParam(1,1) = ratioVals(i)*.5/3600;
%     coupleParam(2,1) = (1-ratioVals(i))*.5/3600;
%     paramList = [9.5582, 3.2924, 1.2876, 0.0499, 0.6165, 0.4856, 11.0546, 3.0365, 1.2880, 0.6255, 0.5021, coupleParam(1,:), coupleParam(2,:), 1];
    paramList = pSol;
    actinInhib =  1;%1 / (1 + (inhibMag(i)/paramList(18)));
    ROCKInhib = 1;% / (1 + (inhibMag(i)/2));
    MRTFInhib = 1;%inhibMag(i);
    YAPInhib = 1;%inhibMag(i);
    cytoDConc = 0;%inhibMag(i);
    LATSFactor = inhibMag(i);%7.5;
    inhibVec = [actinInhib, ROCKInhib, MRTFInhib, YAPInhib, cytoDConc, LATSFactor]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)] 
    [T,Y, ySS] = MechanoCircadianModel([0 maxTime], [stiffnessVals(i),inf,inhibMag(i)], paramList, inhibVec, 0);%kraMult(i));
%     yyaxis left
%     plot(T/3600 - 48,Y(:,5),'-','LineWidth',1)
%     ylabel('PER protein')
%     hold on
%     xlabel('Time [h]')
% %     yyaxis right
%     plot(T/3600 - 48,Y(:,4),'-')
% %     ylabel('BMAL1 protein')
% %     hold on
%     xlim([0 96])
%     prettyGraph
    oscVarIdx = 2;
    if mod(i-1,1)==0
        plot(T/(24*3600) - 1,Y(:,oscVarIdx),'LineWidth',2)%,'Color',colorSeries(i,:))
        xlim([0 5])
        ylabel('PER abundance')
        xlabel('Time (days)')
        prettyGraph
    end
    posTimeIdx = find(T>0);%T>48*3600);
    [pks,locs] = findpeaks(Y(posTimeIdx,oscVarIdx),'MinPeakProminence',.01);
    [troughs,troughLocs] = findpeaks(-Y(posTimeIdx,oscVarIdx),'MinPeakProminence',.01);
    eqTimeIdx = find(T>72*3600);
    YAPTAZEq(i) = ySS(15)/(ySS(17)+ySS(18));%mean(Y(eqTimeIdx,2));
    MRTFEq(i) = ySS(25)/ySS(26);%mean(Y(eqTimeIdx,3));
    FcytoEq(i) = ySS(5);
    GactinEq(i) = ySS(9);
    if ~isempty(pks)
        TShift = T(posTimeIdx);
        period(i) = mean(diff(TShift(locs(2:end))));
        amplitude(i) = mean(pks)-mean(-troughs);
    else
        period(i) = nan;
        amplitude(i) = nan;
    end
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

%% generate population
% lower YAPTAZ in nucleus: 1.25 uM
% higher YAPTAZ in nucleus: 2 uM
fixedParam = [2, 8, 14, 17]; % don't allow exponents to vary
popVar = 0.2;
varVec = popVar*ones(20,1);
varVec(fixedParam) = 0;
numCells = 200;
figure
hold on
maxTime = 3600*1000;
tInterp = 0:960;
oscStored = zeros(numCells,length(tInterp));
YAPTAZStored = zeros(numCells,1);
MRTFStored = zeros(numCells,1);
oscVarIdx = 1; % look at BMAL1 here
for i = 1:numCells
    varVecCur = varVec.*randn(size(varVec));
    pCur = pSol' .* exp(varVecCur(1:19));
    inhibVec = ones(1,4); %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)]
    inhibVec(4) = 1;%0.1; %10% phosphorylation
    inhibVec(5) = 0; % 0 cytD
    inhibVec(6) = 1; %high density LATS factor=7.5, low density=1
    stiffness = 1;
    [T,Y, ySS] = MechanoCircadianModel([0 maxTime], [stiffness,inf,0], pCur, inhibVec, popVar);
    plot(T/3600 - 48, Y(:,oscVarIdx))
    oscStored(i,:) = interp1(T/3600, Y(:,oscVarIdx), tInterp);
    YAPTAZStored(i) = ySS(15);
    MRTFStored(i) = ySS(25);
end

% compute reference period for zero variance
[T,Y, ySS] = MechanoCircadianModel([0 maxTime], [1e5,inf,0], pSol, inhibVec, 0);
posTimeIdx = find(T>0);
[pks,locs] = findpeaks(Y(posTimeIdx,oscVarIdx),'MinPeakProminence',.01);
TShift = T(posTimeIdx);
refPeriod = mean(diff(TShift(locs)));

%% "Circadian power fraction"
% maxTime = 500*3600;
% [T,Y] = MechanoCircadianModel([0 maxTime], [100,inf,0], pSol);%kraMult(i));
tInterp = 48*3600:3600:960*3600;
figure
hold on
Fs = 1/3600;
N = size(oscStored, 2);
allPower = zeros(size(oscStored,1), length(0:Fs/N:Fs/2)); 
for i = 1:size(oscStored,1)
    % Fs = 1; % 1 hr sampling interval
    % tInterp = 0:(1/Fs):maxTime;
    % yInterp = interp1(tOsc,oscStored(i,:),tInterp);
    yInterp = oscStored(i,48:end);
%     hoursSmooth = 48;
%     yWindow = smooth(yInterp, round(hoursSmooth*3600*Fs));
    yWindow = mean(yInterp);
    yInterp = yInterp' - yWindow;
%     % eliminate edge effects
%     keepLogic = tInterp > hoursSmooth*3600 & tInterp < (tInterp(end)-hoursSmooth*3600);
%     tInterp = tInterp(keepLogic);
%     yInterp = yInterp(keepLogic);
    Y_fft = fft(yInterp);
    Y_fft = Y_fft(1:floor(N/2+1));
    curPower = (1/(Fs*N)) * abs(Y_fft).^2;
    curPower(2:end-1) = 2*curPower(2:end-1);
    freq = 0:Fs/N:Fs/2;
    plot(freq*3600, curPower)
    xlim([0 0.4])
    allPower(i,:) = curPower;
end
avgPower = mean(allPower,1);
plot(freq*3600,avgPower,'LineWidth',2)
[~,maxIdx] = max(avgPower);
peakFreq = 1/(refPeriod);%freq(maxIdx);
powerFraction = zeros(size(allPower,1),1);
for i = 1:size(allPower,1)
    totPower = trapz(freq, allPower(i,:));
    freqInterp = peakFreq - 0.2/(3600*24):0.001/(3600*24):peakFreq + 0.2/(3600*24);
    powerInterp = interp1(freq, allPower(i,:), freqInterp);
    powerFraction(i) = trapz(freqInterp, powerInterp) / totPower;
end

%% Fit YAPTAZ circadian model to stiffness data
p0 = [8*3600; 2.5; 1/3600; .04; 0.5/3600; 0.4/3600; 8*3600; 2.5; 1/3600; .5; .4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 1; 3.25];
lowerLim = p0/4;
upperLim = p0*4;
fixedParam = [2, 7, 8, 14, 17];
for i = fixedParam
    lowerLim(i) = p0(i);
    upperLim(i) = p0(i);
end
expParam = [2,8,14,17]; 
for i = expParam
    if lowerLim(i)<1
        lowerLim(i) = 1;
    end
end
options = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon);%, 'MaxTime', 60*20);
pSol = particleswarm(@pToObj_CircadianClock, length(p0), lowerLim, upperLim, options);%p0);
pToObj_CircadianClock(pSol)
save('pSol_justPeriodWrongTrend_var4x.mat','pSol')

%% Fit MRTF
MRTFParamInit = [1.5e6; 1; 10; 1; 100; 0.0461];
lowerLim = MRTFParamInit/10;
upperLim = MRTFParamInit*10;
% lowerLim([1]) = MRTFParamInit([1]);
% upperLim([1]) = MRTFParamInit([1]);
options = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon, 'FunctionTolerance', 1e-9);%, 'MaxTime', 60*20);
pSolMRTF = particleswarm(@pToObj_MRTF, length(MRTFParamInit), lowerLim, upperLim, options);%p0);
% pSolMRTF = lsqnonlin(@pToF_MRTF, MRTFParamInit, lowerLim, upperLim);
pToObj_MRTF(pSolMRTF)

%% plot MRTF cases
stiffnessVals = logspace(-1,4,100);
inhibVals = linspace(0,2,100);
MRTFEq = zeros(size(stiffnessVals));
YAPTAZEq = zeros(size(stiffnessVals));
FactinEq = zeros(size(stiffnessVals));
GactinEq = zeros(size(stiffnessVals));
for i = 1:length(stiffnessVals)
    inhibVec = [1,1,1,1, 0];
    stiffnessVec = [stiffnessVals(i),inf,0];
%     stiffnessVec = [1e5,inf,0];
    SSVar = MechanoSS(stiffnessVec, inhibVec, pSolMRTF);
    MRTFEq(i) = SSVar(25)/SSVar(26);
    YAPTAZEq(i) = SSVar(15)/(SSVar(17)+SSVar(18));
    FactinEq(i) = SSVar(5);
    GactinEq(i) = SSVar(9);
end
% semilogx(stiffnessVals,YAPTAZEq)
% hold on
% semilogx(stiffnessVals,MRTFEq)
% hold on

function obj = pToObj_CircadianClock(p)
    maxTime = 3600*1000;
%     latruncTreatments = [0, 1, 2];%, 5];
    cytoDTreatments = [0,1,2];
    stiffnessTests = [1e4, 300, 19];
%     periodVec = 3600*[24.9, 25.5, 26];%, 24.6];
%     amplVec = [1, 2/1.2, 3/1.2];%, 3.3/1.2];
%     periodVec = 3600*[24, 25, 25];%, 24.6];
%     amplVec = [1, 1.6/1.2, 2/1.2];%, 2.3/1.2];
    periodVec1 = 3600*[23.8, 24.6, 25.6];
    periodVec2 = 3600*[23.8, 24, 24.8];
    periodVec = (periodVec1 + periodVec2)/2;
    amplVec1 = [1, 105/45, 85/45];
    amplVec2 = [1, 145/80, 95/80];
    amplVec = (amplVec1 + amplVec2) / 2;
    periodTest = zeros(size(periodVec1));
    amplTest = zeros(size(amplVec1));
    for i = 1:length(stiffnessTests)
%         actinPolymBlock =  1 / (1 + (latruncTreatments(i)/p(18)));
        actinPolymBlock = 1;
        inhibVec = [actinPolymBlock, 1, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)] 
        [t,y] = MechanoCircadianModel([0 maxTime], [stiffnessTests(i),inf,0], p, inhibVec, 0);
        if any(~isreal(y))
            obj = 1e6;
            return
        end
        [pks,locs] = findpeaks(y(:,2),'MinPeakProminence',.01);
        [troughs,~] = findpeaks(-y(:,2),'MinPeakProminence',.01);
        if length(pks)>10 && length(troughs)>10
            lastExtremum = min([length(pks),length(troughs)]); 
            amplitudeStored = pks(2:lastExtremum) - (-troughs(2:lastExtremum));
            amplTest(i) = mean(amplitudeStored);
            if amplitudeStored(end) < 0.8*amplitudeStored(1) % cannot be decaying rapidly (not sustained osc)
                periodTest(i) = t(end);
            else
                periodTest(i) = mean(diff(t(locs(2:end))));
            end
        else
            periodTest(i) = t(end);
            amplTest(i) = mean(y(:,2));
        end
    end
    amplTest = amplTest/amplTest(1); % normalize to control amplitude
    obj = sum((periodTest-periodVec1(end:-1:1)).^2); % just period, reversed
%     obj = sum(50*((periodTest - periodVec1)./periodVec1).^2 + ((amplTest - amplVec1)./amplVec1).^2);% +...
%         sum(50*((periodTest - periodVec2)./periodVec2).^2 + ((amplTest - amplVec2)./amplVec2).^2);
%     obj = sum(50*((periodTest - periodVec)./periodVec).^2 + ((amplTest - amplVec)./amplVec).^2);
end

function obj = pToObj_MRTF(p)
    stiffnessTests = [1, 1, 1e5, 1e5];
    cytDTests = [0, 2.5, 0, 2.5];
    MRTFEqData = [0.5, 0.8, 1, 4];
    MRTFEqTest = zeros(size(MRTFEqData));
    for i = 1:length(cytDTests)
        inhibVec = [1,1,1,1, cytDTests(i)];
        stiffnessVec = [stiffnessTests(i),inf,0];
        SSVar = MechanoSS(stiffnessVec, inhibVec, p);
        MRTFEqTest(i) = SSVar(25)/SSVar(26);
    end
    obj = sum((MRTFEqData-MRTFEqTest).^2);
end

function F = pToF_MRTF(p)
    stiffnessTests = [1, 1, 1e5, 1e5];
    cytDTests = [0, 2.5, 0, 2.5];
    MRTFEqData = [0.5, 0.8, 1, 4];
    MRTFEqTest = zeros(size(MRTFEqData));
    for i = 1:length(cytDTests)
        inhibVec = [1,1,1,1, cytDTests(i)];
        stiffnessVec = [stiffnessTests(i),inf,0];
        SSVar = MechanoSS(stiffnessVec, inhibVec, p);
        MRTFEqTest(i) = SSVar(25)/SSVar(26);
    end
    F = abs(MRTFEqData - MRTFEqTest);
end