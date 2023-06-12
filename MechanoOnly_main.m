%% run sims with diff stiffnesses
stiffnessVals = logspace(0,2,41);
maxTime = 3600*120;
tauVals = zeros([24,length(stiffnessVals)]);
figure
hold on
for i = 1:length(stiffnessVals)
    [T,Y,TRed,SSVarMat,tauValsCur] = MechanoOnlyModel([0 maxTime],[stiffnessVals(i),3600*240]);
    tauVals(:,i) = tauValsCur;
    if ~mod(i-1,10)
        plot(T/3600,Y(:,15)./(Y(:,17)+Y(:,18)));
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