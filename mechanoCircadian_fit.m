%% Fit YAPTAZ circadian model to stiffness data
p0 = [12*3600; 2; 1/3600; .04; 1/3600; 0.4/3600; 6*3600; 2; 1/3600; .5; .4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; 0.4/3600];
% p0(end) = 10;
lowerLim = p0/4;
upperLim = p0*4;
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
options = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,'PlotFcn','pswplotbestf',...
    'FunctionTolerance',1e-4,'PlotFcn','pswplotbestf','MaxStallIterations',15);%, 'MaxTime', 60*20);
pSol = particleswarm(@pToObj_CircadianClockFibroblast, length(p0), lowerLim, upperLim, options);%p0);
pToObj_CircadianClockFibroblast(pSol)
save('pSol_var4x_4xDynand1xPeriod_lucModel_2delays.mat','pSol')

%%
pToObj_CircadianClockFibroblast(pSol)

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
figure
sp3 = subplot(1,3,3);
xlabel('Time (days)')
ylabel('PER/CRY expression rate (nM/hr)')
prettyGraph
hold on
stiffnessTests = logspace(0, 7.5, 50);
expStiffness = [1e7, 300, 19];
stiffnessLoop = [stiffnessTests, expStiffness];
periodStiffness = zeros(size(stiffnessTests));
amplStiffness = zeros(size(stiffnessTests));
for i = 1:length(stiffnessLoop)
    inhibVec = [1, 1, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [curPeriod, curAmpl, t, y, rawOutput] = conditionToOutputs(pSol, stiffnessLoop(i), inhibVec, 3600*480);
    if i > length(stiffnessTests)
        oscDynamics = y(:,3);%3600*pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
        plot(t/(24*3600),oscDynamics,'LineWidth',1)
        xlim([0 5])
%         ylim([0 .25])
    else
        periodStiffness(i) = curPeriod;
        amplStiffness(i) = curAmpl;
    end
end
legend('Control (glass)','300 kPa PEKK scaffold','19 kPa PCL scaffold')
sp1 = subplot(1,3,1);
semilogx(stiffnessTests, periodStiffness/3600, 'Color', 'b', 'LineWidth', 1)
hold on
errorbar(expStiffness, FibroblastStiffnessPeriod(:,1),...
    FibroblastStiffnessPeriod(:,2), FibroblastStiffnessPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('Stiffness (kPa)')
ylabel('Period (hr)')
legend('Model','Experiments')
ylim([20 28])
prettyGraph
sp2 = subplot(1,3,2);
semilogx(stiffnessTests, amplStiffness/amplStiffness(end), 'Color', 'b', 'LineWidth', 1)
hold on
errorbar(expStiffness, FibroblastStiffnessAmpl(:,1),...
    FibroblastStiffnessAmpl(:,2), FibroblastStiffnessAmpl(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('Stiffness (kPa)')
ylabel('Normalized amplitude')
legend('Model','Experiments')
ylim([0 2.5])
prettyGraph
pos1 = get(sp1, 'Position');
pos2 = get(sp2, 'Position');
pos3 = get(sp3, 'Position');
set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])

%% Compare ROCK data
figure
ROCKInhibTests = 0:30;
FibroblastROCKInhib = [0, 10, 20];
ROCKInhibLoop = [ROCKInhibTests, FibroblastROCKInhib];
periodROCKInhib = zeros(size(ROCKInhibTests));
amplROCKInhib = zeros(size(ROCKInhibTests));
sp3 = subplot(1,3,3);
xlabel('Time (days)')
ylabel('PER/CRY expression rate (nM/hr)')
prettyGraph
hold on
for i = 1:length(ROCKInhibLoop)
    ROCKLevel =  1 / (1 + (ROCKInhibLoop(i)/pSol(18)));
    inhibVec = [1, ROCKLevel, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [curPeriod, curAmpl, t, y] = conditionToOutputs(pSol, 1e7, inhibVec, 480*3600);
    if i > length(ROCKInhibTests)
        oscDynamics = 3600*pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
        plot(t/(24*3600),oscDynamics,'LineWidth',1)
        xlim([0 5])
        ylim([0 .25])
    else
        periodROCKInhib(i) = curPeriod;
        amplROCKInhib(i) = curAmpl;
    end
end
legend('Control','10uM Y27632','20uM Y27632')
sp1 = subplot(1,3,1);
plot(ROCKInhibTests, periodROCKInhib/3600,'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastROCKInhib, FibroblastROCKInhibPeriod(:,1),...
    FibroblastROCKInhibPeriod(:,2), FibroblastROCKInhibPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlim([-1,25])
ylim([20 28])
xlabel('Y27632 concentration (uM)')
ylabel('Period (hr)')
prettyGraph
sp2 = subplot(1,3,2);
plot(ROCKInhibTests, amplROCKInhib/amplROCKInhib(1),'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastROCKInhib, FibroblastROCKInhibAmpl(:,1),...
    FibroblastROCKInhibAmpl(:,2), FibroblastROCKInhibAmpl(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('Y27632 concentration (uM)')
ylabel('Normalized amplitude')
xlim([-1,25])
ylim([0 2.5])
prettyGraph
pos1 = get(sp1, 'Position');
pos2 = get(sp2, 'Position');
pos3 = get(sp3, 'Position');
set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])

%% CytD Tests
figure
sp3 = subplot(1,3,3);
xlabel('Time (days)')
ylabel('PER/CRY expression rate (nM/hr)')
prettyGraph
hold on
CytDTests = 0:.2:6;
FibroblastCytD = [0, 1, 2, 5];
CytDLoop = [CytDTests, FibroblastCytD];
periodCytD = zeros(size(CytDTests));
amplCytD = zeros(size(CytDTests));
for i = 1:length(CytDLoop)
    inhibVec = [1, 1, 1, 1, CytDLoop(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [curPeriod, curAmpl, t, y] = conditionToOutputs(pSol, 1e7, inhibVec, 480*3600);
    if i > length(CytDTests)
        oscDynamics = 3600*pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
        plot(t/(24*3600),oscDynamics,'LineWidth',1)
        xlim([0 5])
        ylim([0 .25])
    else
        periodCytD(i) = curPeriod;
        amplCytD(i) = curAmpl;
    end
end
legend('Control','1uM CytD','2uM CytD','5uM CytD')
sp1 = subplot(1,3,1);
plot(CytDTests, periodCytD/3600,'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastCytD, FibroblastCytDPeriod(:,1),...
    FibroblastCytDPeriod(:,2), FibroblastCytDPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('CytD concentration (uM)')
ylabel('Period (hr)')
xlim([-.2 6])
ylim([20 28])
prettyGraph
sp2 = subplot(1,3,2);
plot(CytDTests, amplCytD/amplCytD(1),'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastCytD, FibroblastCytDAmpl(:,1),...
    FibroblastCytDAmpl(:,2), FibroblastCytDAmpl(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('CytD concentration (uM)')
ylabel('Normalized amplitude')
xlim([-.2 6])
ylim([0 2.5])
prettyGraph
pos1 = get(sp1, 'Position');
pos2 = get(sp2, 'Position');
pos3 = get(sp3, 'Position');
set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])

%% LatB Tests
figure
sp3 = subplot(1,3,3);
xlabel('Time (days)')
ylabel('PER/CRY expression rate (nM/hr)')
prettyGraph
hold on
LatBTests = 0:.1:2.5;
FibroblastLatB = [0, 1, 2];
LatBLoop = [LatBTests, FibroblastLatB];
periodLatB = zeros(size(LatBTests));
amplLatB = zeros(size(LatBTests));
for i = 1:length(LatBLoop)
    kra_cur = 1 / (1 + (LatBLoop(i)/pSol(26)));
    inhibVec = [kra_cur, 1, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [curPeriod, curAmpl, t, y] = conditionToOutputs(pSol, 1e7, inhibVec, 480*3600);
    if i > length(LatBTests)
        oscDynamics = 3600*pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
        plot(t/(24*3600),oscDynamics,'LineWidth',1)
        xlim([0 5])
        ylim([0 1])
    else
        periodLatB(i) = curPeriod;
        amplLatB(i) = curAmpl;
    end
end
legend('Control','1uM LatB','2uM LatB')
sp1 = subplot(1,3,1);
plot(LatBTests, periodLatB/3600,'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastLatB, FibroblastLatBPeriod(:,1),...
    FibroblastLatBPeriod(:,2), FibroblastLatBPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('LatB concentration (uM)')
ylabel('Period (hr)')
prettyGraph
xlim([-.1 2.5])
ylim([20 28])
sp2 = subplot(1,3,2);
plot(LatBTests, amplLatB/amplLatB(1),'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastLatB, FibroblastLatBAmpl(:,1),...
    FibroblastLatBAmpl(:,2), FibroblastLatBAmpl(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('LatB concentration (uM)')
ylabel('Normalized amplitude')
xlim([-.1 2.5])
ylim([0 2.5])
prettyGraph
pos1 = get(sp1, 'Position');
pos2 = get(sp2, 'Position');
pos3 = get(sp3, 'Position');
set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])

%% Jas Tests
figure
sp3 = subplot(1,3,3);
xlabel('Time (days)')
ylabel('PER/CRY expression rate (nM/hr)')
prettyGraph
hold on
JasTests = 0:.02:.6;
FibroblastJas = [0, 0.1, 0.2, 0.5];
JasLoop = [JasTests,FibroblastJas];
periodJas = zeros(size(JasTests));
amplJas = zeros(size(JasTests));
for i = 1:length(JasLoop)
    kra_cur = 1 + pSol(27)*JasLoop(i) / (pSol(28) + JasLoop(i));
    inhibVec = [kra_cur, 1, 1, 1, 0]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytoDConc] 
    [curPeriod, curAmpl, t, y] = conditionToOutputs(pSol, 1e7, inhibVec);
     if i > length(JasTests)
        oscDynamics = 3600*pSol(9)./(1 + (pSol(10)./y(:,1)).^pSol(8));
        plot(t/(24*3600),oscDynamics,'LineWidth',1)
        xlim([0 5])
%         ylim([.1 .3])
    else
        periodJas(i) = curPeriod;
        amplJas(i) = curAmpl;
    end
end
legend('Control','0.1uM Jas','0.2uM Jas','0.5uM Jas')
sp1 = subplot(1,3,1);
plot(JasTests, periodJas/3600,'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastJas, FibroblastJasPeriod(:,1),...
    FibroblastJasPeriod(:,2), FibroblastJasPeriod(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('Jas concentration (uM)')
ylabel('Period (hr)')
xlim([-.02 .6])
ylim([20 28])
prettyGraph
sp2 = subplot(1,3,2);
plot(JasTests, amplJas/amplJas(1),'LineWidth',1,'Color','b')
hold on
errorbar(FibroblastJas, FibroblastJasAmpl(:,1),...
    FibroblastJasAmpl(:,2), FibroblastJasAmpl(:,2),...
    'LineStyle','none','LineWidth',1,'Marker','s','Color','r')
xlabel('Jas concentration (uM)')
ylabel('Normalized amplitude')
xlim([-.02 .6])
ylim([0 2.5])
prettyGraph
pos1 = get(sp1, 'Position');
pos2 = get(sp2, 'Position');
pos3 = get(sp3, 'Position');
set(sp1, 'Position', [pos1(1), pos1(2)+.1*pos1(4), 0.8*pos1(3), 0.9*pos1(4)])
set(sp2, 'Position', [pos2(1)-0.2*pos1(3), pos2(2)+.1*pos2(4), 0.8*pos2(3), 0.9*pos2(4)])
set(sp3, 'Position', [pos3(1)-0.4*pos2(3), pos3(2)+.1*pos3(4), pos3(3)*1.5, 0.9*pos3(4)])

%% assess goodness of fit w.r.t. period
periodVec = [23.8447, 24.5874, 25.5922,...
             24.7122, 25.1024, 25.1707,...
             24.2103, 24.9957, 24.9700, 24.6352,...
             24.9389, 25.4105, 25.9694,...
             23.6788, 23.2228, 22.3938, 20.4041];
[~,periodTest] = pToObj_CircadianClockFibroblast(pSol);
gof = goodnessOfFit(periodTest', periodVec', 'MSE')

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
        kra_cur = 1 / (1 + (LatBTests(i)/p(26))) + p(27)*JasTests(i) / (p(28) + JasTests(i));
        inhibVec = [kra_cur, ROCKLevel, 1, 1, CytDTests(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC), CytD] 
        [periodTest(i), amplTest(i), t, y] = conditionToOutputs(p, stiffnessTests(i), inhibVec);
        oscDynamics = y(:,3);%3600*p(9)./(1 + (p(10)./y(:,1)).^p(8));
        modelStored{i} = [t, oscDynamics];
        minVals(i) = min(oscDynamics);
    end
    periodTest([4,7,11,14]) = periodTest(1); % all control cases
    periodTest = periodTest/3600;
    amplTestRaw = amplTest;
    amplTest = amplTest/amplTest(1); % normalize to control amplitude
    amplTest([4,7,11,14]) = amplTest(1); % all control cases
    modelStored{4} = modelStored{1};
    modelStored{7} = modelStored{1};
    modelStored{11} = modelStored{1};
    modelStored{14} = modelStored{1};

    obj = 50*(sum(4*((periodTest(1:3) - periodVec(1:3))./periodVec(1:3)).^2) + ...
              sum(3*((periodTest(4:6) - periodVec(4:6))./periodVec(4:6)).^2) + ...
              sum(3*((periodTest(7:10) - periodVec(7:10))./periodVec(7:10)).^2) + ...
              sum(3*((periodTest(11:13) - periodVec(11:13))./periodVec(11:13)).^2) + ...
              sum(3*((periodTest(14:17) - periodVec(14:17))./periodVec(14:17)).^2));% +...
%           .25*(sum(4*((amplTest(1:3) - amplVec(1:3))./amplVec(1:3)).^2) + ...
%              sum(3*((amplTest(4:6) - amplVec(4:6))./amplVec(4:6)).^2) + ...
%              sum(3*((amplTest(7:10) - amplVec(7:10))./amplVec(7:10)).^2) + ...
%              sum(3*((amplTest(11:13) - amplVec(11:13))./amplVec(11:13)).^2) + ...
%              sum(3*((amplTest(14:17) - amplVec(14:17))./amplVec(14:17)).^2));
    nVec = [4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
%     figure
    for i = 1:length(periodVec)
        yNormalized = (modelStored{i}(:,2) - min(minVals))/amplTestRaw(1);
        tTest = modelStored{i}(:,1);
        PERTest = 0.5*amplVec(i)*(cos(2*pi*tTest/(periodVec(i)*3600))) + max(amplVec);
        obj = obj + nVec(i)*200e-4*sum(((yNormalized - PERTest)./amplVec(i)).^2);
%         subplot(2,9,i)
%         title(sprintf("%d", i))
%         plot(tTest, yNormalized)
%         hold on
%         plot(tTest, PERTest)
    end
end
