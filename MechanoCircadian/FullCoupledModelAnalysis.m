%% quick phase diagram from full coupled Circadian-mechanotransduction model
coupingParam = zeros(1,2);
couplingParam(1) = (1.0 * 0.001 / 3600.0); % nu_YTB
couplingParam(2) = (1.0 / (1000.0 * 3600.0)); % nu_MRTF
stiffnessVal = 1e5; % TCP
nu_YTBVals = (0:.1:4)*couplingParam(1);
nu_MRTFVals = (0:.1:4)*couplingParam(2);
[nu_YTBMesh, nu_MRTFMesh] = meshgrid(nu_YTBVals,nu_MRTFVals);
periodMesh = zeros(size(nu_YTBMesh));
amplMesh = zeros(size(nu_YTBMesh));
mPER_avg = zeros(size(nu_YTBMesh));
mCRY_avg = zeros(size(nu_YTBMesh));
mBMAL1_avg = zeros(size(nu_YTBMesh));
tRange = [0, 3600*500];
oscVarIdx = 42;
for i = 1:length(nu_YTBMesh(:))
    couplingCur = zeros(2,4);
    couplingCur(1,1) = nu_YTBMesh(i);
    couplingCur(2,2:3) = nu_MRTFMesh(i);
    [T,Y] = YAPTAZ_FullClockNew(tRange,[],[],stiffnessVal,couplingCur);
    [pks,locs] = findpeaks(Y(:,oscVarIdx),'MinPeakProminence',.01*1e-3);
    [troughs,troughLocs] = findpeaks(-Y(:,oscVarIdx),'MinPeakProminence',.01*1e-3);
    mPER_avg(i) = trapz(T,Y(:,17))/range(T);
    mCRY_avg(i) = trapz(T,Y(:,39))/range(T);
    mBMAL1_avg(i) = trapz(T,Y(:,12))/range(T);
    if length(pks)>4
        periodMesh(i) = mean(diff(T(locs)));
        amplMesh(i) = mean(pks)-mean(-troughs);
    else
        periodMesh(i) = nan;
        amplMesh(i) = nan;
    end
end

%% Fit values for desired shift in periods
% couplingRef = 1.0 / (1000.0 * 3600.0);
% couplingParam0 = 0.5*couplingRef*ones(8,1);
lowerLim = 0.8*ones(108,1);
upperLim = 1.25*ones(108,1);
lowerLim(1:8) = 0;
upperLim(1:8) = 5;
% could restrict coupling to postulated connections from papers
upperLim(2:3) = 0; % no YAP/TAZ coupling to PER/CRY
upperLim([5,8]) = 0; % no MRTF coupling to BMAL or REVERB
options = optimoptions('particleswarm','UseParallel', false,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','FunctionTolerance',1e-3,'MaxStallIterations',10);%, 'MaxTime', 60*20);
pSol = particleswarm(@pToObj_FullCoupling, length(lowerLim), lowerLim, upperLim, options);%p0);
pToObj_FullCoupling(pSol)

%% plot different stiffnesses
circad_param = pSol;
% couplingSol = reshape(pSol,[4,2])';
stiffnessVals = logspace(-1,3,20);
tRange = [0, 3600*500];
oscVarIdx = 18;
period = zeros(size(stiffnessVals));
ampl = zeros(size(stiffnessVals));
% mPER_avg = zeros(size(stiffnessVals));
% mCRY_avg = zeros(size(stiffnessVals));
% mBMAL1_avg = zeros(size(stiffnessVals));
figure
hold on
for i = 1:length(stiffnessVals)
    inhibVec = [1,1,1,1,0];
    stiffnessVec = [stiffnessVals(i),inf,0];
    MRTFParam = [1, 3.82752594304687, 20.9581804455536, 0.366250664261199, 86.6114163570272, 0.184399999956306];
    SSVar = MechanoSS(stiffnessVec, inhibVec, MRTFParam);
    [T,Y] = CircadianOnly_full(tRange, circad_param, SSVar([15,26]));
    [pks,locs] = findpeaks(Y(:,oscVarIdx),'MinPeakProminence',.01*1e-3);
    [troughs,troughLocs] = findpeaks(-Y(:,oscVarIdx),'MinPeakProminence',.01*1e-3);
    % mPER_avg(i) = trapz(T,Y(:,17))/range(T);
    % mCRY_avg(i) = trapz(T,Y(:,39))/range(T);
    % mBMAL1_avg(i) = trapz(T,Y(:,12))/range(T);
    if length(pks)>4
        period(i) = mean(diff(T(locs)));
        ampl(i) = mean(pks)-mean(-troughs);
    else
        period(i) = nan;
        ampl(i) = nan;
    end
    plot(T/3600, Y(:,oscVarIdx))
end

function obj = pToObj_FullCoupling(circad_param)
    maxTime = 3600*500;
    stiffnessTests = [1e4, 300, 19];
    periodVec1 = 3600*[23.8, 24.6, 25.6];
    periodVec2 = 3600*[23.8, 24, 24.8];
    periodVec = (periodVec1 + periodVec2)/2;
    amplVec1 = [1, 105/45, 85/45];
    amplVec2 = [1, 145/80, 95/80];
    amplVec = (amplVec1 + amplVec2) / 2;
    periodTest = zeros(size(periodVec1));
    amplTest = zeros(size(amplVec1));
    for i = 1:length(stiffnessTests)
        inhibVec = [1,1,1,1,0];
        stiffnessVec = [stiffnessTests(i),inf,0];
        MRTFParam = [1, 3.82752594304687, 20.9581804455536, 0.366250664261199, 86.6114163570272, 0.184399999956306];
        SSVar = MechanoSS(stiffnessVec, inhibVec, MRTFParam);
        [t,y] = CircadianOnly_full([0 maxTime], circad_param, SSVar([15,26]));
        if any(~isreal(y))
            obj = 1e6;
            return
        end
        [pks,locs] = findpeaks(y(:,18),'MinPeakProminence',.01*1e-3);
        [troughs,~] = findpeaks(-y(:,18),'MinPeakProminence',.01*1e-3);
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
    obj = sum(50*((periodTest - periodVec1)./periodVec1).^2 + ((amplTest - amplVec1)./amplVec1).^2) +...
        sum(50*((periodTest - periodVec2)./periodVec2).^2 + ((amplTest - amplVec2)./amplVec2).^2);
%     obj = sum(50*((periodTest - periodVec)./periodVec).^2 + ((amplTest - amplVec)./amplVec).^2);
end