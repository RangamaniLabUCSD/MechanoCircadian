%% Plot mechanotransduction dynamics
maxTime = 3600*6;
[T,Y] = MechanoOnlyModel([0 maxTime],[100,inf]);
figure
% subplot(1,2,1)
hold on
plot(T/60,Y(:,15)./(Y(:,17)+Y(:,18)),'LineWidth',1);
plot(T/60,Y(:,26)./(Y(:,25)+Y(:,26)),'LineWidth',1);
xlim([0 2*60])
legend('YAP/TAZ','MRTF')
ylabel('Nuclear to cytoplasmic ratio')
xlabel('Time (minutes)')
prettyGraph
% subplot(1,2,2)
% hold on
% plot(T/60,Y(:,5)/Y(1,5),'LineWidth',1)
% plot(T/60,Y(:,21)/Y(1,21),'LineWidth',1)
% xlim([0 2*60])
% legend('F-actin','Active myosin')
% prettyGraph

%% Test reduced YAP-TAZ system dynamics
stiffnessVals = logspace(0,2,41);
maxTime = 3600*120;
tauVals = zeros([24,length(stiffnessVals)]);
figure
hold on
for i = 1:length(stiffnessVals)
    [T,Y,TRed,SSVarMat,tauValsCur] = MechanoOnlyModel([0 maxTime],[stiffnessVals(i),inf]);
    tauVals(:,i) = tauValsCur;
    if ~mod(i-1,10)
        plot(T/60,Y(:,15)./(Y(:,17)+Y(:,18)));
        plot(TRed/3600,SSVarMat(:,15)./(SSVarMat(:,17)+SSVarMat(:,18)),'--')
%         plot(TRed/3600,(SSVals(15)/(SSVals(17)+SSVals(18)))*ones(size(TRed)),'--')
    end
end

figure
varsIdx = [16,6,7,11,21,24,13,5,4,8,17,15];
namesList = {'Fakp','RhoAGTP','mDiaA','ROCKA','MyoA','LIMKA',...
                'CofilinNP','Fcyto','LaminA','NPCA','YAPTAZP','YAPTAZnuc'};
for i = varsIdx
    loglog(stiffnessVals,tauVals(i,:))
    hold on
end
legend(namesList)

%% phase plane analysis
stiffnessVal = 100;
maxTime = 1e4;
[T,Y,TRed,SSVarMat,tauValsCur,param] = MechanoOnlyModel([0 maxTime],[stiffnessVal,inf]);

saveVid = true;
if saveVid
    frameRate = 25;
    path = uigetdir('Choose where to save video');
    myVideo = VideoWriter(fullfile(path,'100kPaPhasePlane.mp4'),'MPEG-4');
    myVideo.FrameRate = frameRate;
    open(myVideo)
end

%phase plane
YPhase = (0:.4:4)';
LPhase = (0:100:3500)';
[LGrid,YGrid] = meshgrid(LPhase,YPhase);
YPhase = YGrid(:);
LPhase = LGrid(:);
figure
for t = 0:10:5000
    clf
    dydt = zeros([length(YPhase),2]);
    for i = 1:length(LPhase)
        dydt(i,:) = YAPTAZPhasePlane(t,[LPhase(i),YPhase(i)],param);
    end
    dydt(:,1) = dydt(:,1)/3500;
    dydt(:,2) = dydt(:,2)/4;
    arrowMag = sqrt(dydt(:,1).^2 + dydt(:,2).^2);
    quiver(LPhase/3500,YPhase/4,dydt(:,1),dydt(:,2),'AutoScaleFactor',1)
    hold on
    [~,nc1,nc2,fp] = YAPTAZPhasePlane(t,[LPhase(i),YPhase(i)],param);
    plot(nc1(1,:)/3500,nc1(2,:)/4,'-r')
    plot(nc2(1,:)/3500,nc2(2,:)/4,'-r')
    plot(fp(1)/3500,fp(2)/4,'pr','MarkerSize',8)

    LCur = interp1(TRed,SSVarMat(:,4),t);
    YCur = interp1(TRed,SSVarMat(:,15),t);
    plot(LCur/3500,YCur/4,'Marker','o','MarkerSize',8,'Color','g','LineWidth',1)
    title(sprintf('t = %.1f s',t))
    xlim([0 1])
    ylim([0 1])
    xlabel('Rel LaminA conc')
    ylabel('Rel YAP/TAZ nuclear conc')
    drawnow
    if saveVid
        img = print('-RGBImage');
        %                 img = print('-RGBImage','-r80');
        writeVideo(myVideo, img)
    end
end
if saveVid
    close(myVideo)
end

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
figure
semilogx(stiffnessVals,YAPTAZEq)
hold on
semilogx(stiffnessVals,MRTFEq)
hold on

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