%% two variable DDE model of Circadian oscillations
% in this file, we model the circadian system with baseline shifts in
% expression. Together with ddeBifCircadian.m and the
% mechanotransduction-related functions, this allows us to plot oscillation
% period and amplitude over the YAP/TAZ-MRTF phase plane.

%% first define pSol
p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
if ~exist('myBayesianAnalysis','var')
    error('Load in myBayesianAnalysis first')
end
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP','burnIn', 500)
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
fixedParam = [5, 14, 17, 22, 25, 44, 45];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
pSol = p0;
pSol(varyLogic) = pSol(varyLogic) .* modeVals';
p = [pSol(1:11); 0; 0; pSol(35:38); 0; pSol(43)]; % init KeB2, KeP2, KeR2 to zero

% convert to units of hours
p([1,7,17]) = p([1,7,17])/3600;
p([3,5,6,9,11,14,16]) = p([3,5,6,9,11,14,16])*3600;
pRef = p;

%% here, we test a range of KeB2 and KeP2 associated with different values of YAP/TAZ and MRTF (see MechanoCircadian.m for full model)
nCoupleTest = true; % test different values of the Hill coefficient for mechanotransduction-circadian coupling
KdTest = false; % test different values for decay rates in circadian model
timeSpan = [0 960];
YVals = logspace(-1,1.3,20); %40
MVals = logspace(-1,1.3,20); %40
[YMat, MMat] = meshgrid(YVals, MVals);
nCouple = 2;

CytoConv = 1.3851e6;
NucConv = 3.3122e5;
Ytot = 1.4784e6;
Mtot = 1e6;
YConcMat = YMat*Ytot./(CytoConv + YMat*NucConv);
MConcMat = MMat*Mtot./(CytoConv + MMat*NucConv);
KeB2Mat = 3600*(pSol(12)*YConcMat.^nCouple./(pSol(13)^nCouple+YConcMat.^nCouple) +...
                pSol(23)*MConcMat.^nCouple./(pSol(24)^nCouple+MConcMat.^nCouple));
KeP2Mat = 3600*(pSol(15)*MConcMat.^nCouple./(pSol(16)^nCouple+MConcMat.^nCouple) +...
                pSol(20)*YConcMat.^nCouple./(pSol(21)^nCouple+YConcMat.^nCouple)); 
KeR2Mat = 3600*(pSol(41)*MConcMat.^nCouple./(pSol(42)^nCouple+MConcMat.^nCouple) +...
                pSol(39)*YConcMat.^nCouple./(pSol(40)^nCouple+YConcMat.^nCouple));

if nCoupleTest
    nCoupleVals = [1,1.5,2,3,4.5];
    KdBVals = p(6)*ones(size(nCoupleVals));%[1,0.1,3,10,1,1,1,1,1,1]; % KdB (p(6))
    KdPVals = p(11)*ones(size(nCoupleVals));%[1,1,1,1,.2,.5,2,1,1,1]; % KdP (p(11))
    KdRVals = p(16)*ones(size(nCoupleVals));%[1,1,1,1,1,1,1,.1,10,30]; % KdR (p(16))
elseif KdTest %#ok<UNRCH>
    KdBVals = p(6)*[1,0.1,3,10,1,1,1,1,1,1]; % KdB (p(6))
    KdPVals = p(11)*[1,1,1,1,.2,.5,2,1,1,1]; % KdP (p(11))
    KdRVals = p(16)*[1,1,1,1,1,1,1,.1,10,30]; % KdR (p(16))
    nCoupleVals = 2*ones(size(KdBVals));
end
circStoreCell = cell(size(KdBVals));
for k = 1:length(KdBVals)
    p(6) = KdBVals(k); p(11) = KdPVals(k); p(16) = KdRVals(k);
    nCouple = nCoupleVals(k);
    KeB2Mat = 3600*(pSol(12)*YConcMat.^nCouple./(pSol(13)^nCouple+YConcMat.^nCouple) +...
                pSol(23)*MConcMat.^nCouple./(pSol(24)^nCouple+MConcMat.^nCouple));
    KeP2Mat = 3600*(pSol(15)*MConcMat.^nCouple./(pSol(16)^nCouple+MConcMat.^nCouple) +...
                    pSol(20)*YConcMat.^nCouple./(pSol(21)^nCouple+YConcMat.^nCouple)); 
    KeR2Mat = 3600*(pSol(41)*MConcMat.^nCouple./(pSol(42)^nCouple+MConcMat.^nCouple) +...
                    pSol(39)*YConcMat.^nCouple./(pSol(40)^nCouple+YConcMat.^nCouple));
    period = zeros(size(KeB2Mat));
    oscDecayRate = zeros(size(KeB2Mat));
    amplitude = zeros(size(KeB2Mat));
    amplitudeRatio = zeros(size(KeB2Mat));
    phaseLag = zeros(size(KeB2Mat));
    [numRows, numCols] = size(KeB2Mat);
    for i = 1:numRows
        fprintf('Test %d, row %d of %d\n', k, i, numRows)
        for j = 1:numCols
            pCur = p;
            % figure(timeFig)
            pCur(12) = KeB2Mat(i,j);
            pCur(13) = KeP2Mat(i,j);
            pCur(18) = KeR2Mat(i,j);
            [t2,y2] = DDESolve2(timeSpan,pCur);
            if any(isnan(y2(:,2)))
                y2(:,2) = 0;
            end
            if any(~isreal(y2(:,2)))
                y2(:,2) = real(y2(:,2));
            end
            expressionFunc = pCur(9)./(1+(pCur(10)./y2(:,1)).^pCur(8));
            [period(i,j), amplitude(i,j), oscDecayRate(i,j)] = circOscAnalysis(t2*3600, y2(:,2));
        end
    end
    circStoreCell{k} = {period, amplitude, oscDecayRate};
end

%% plot all phase diagrams after running the KdTest above
figure
YVals = logspace(-1,1.3,20); %40
MVals = logspace(-1,1.3,20); %40
[YMat, MMat] = meshgrid(YVals, MVals);
for i = 2:length(circStoreCell)
    subplot(3,3,i-1)
    oscPeriod = circStoreCell{i}{1}/3600;
    oscDecayRate = circStoreCell{i}{3};
    oscPeriod(oscDecayRate<-.1) = nan; % oscillations decaying too fast to quantify
    surf(YMat, MMat, oscPeriod, 'LineStyle','none','FaceColor','interp')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    daspect([1 1 1])
    view([0 90])
    drawnow
    colorbar
    xlim([.1 20])
    ylim([.1 20])
    clim([min(oscPeriod(oscPeriod>0)) max(oscPeriod(:))])
end

%% compute KeB2 and KeP2 for YAP/TAZ branch
curBranch = BranchesStored{8};
YVals = curBranch{1};
MVals = curBranch{2};
CytoConv = 1.3851e6;
NucConv = 3.3122e5;
Ytot = 1.4784e6;
Mtot = 1e6;
YConcVals = YVals*Ytot./(CytoConv + YVals*NucConv);
MConcVals = MVals*Mtot./(CytoConv + MVals*NucConv);
KeB2Vals = 3600*(pSol(12)*YConcVals.^2./(pSol(13)^2+YConcVals.^2) +...
                pSol(23)*MConcVals.^2./(pSol(24)^2+MConcVals.^2));
KeP2Vals = 3600*(pSol(15)*MConcVals.^2./(pSol(16)^2+MConcVals.^2) +...
                pSol(20)*YConcVals.^2./(pSol(21)^2+YConcVals.^2));

% also compute total baseline expression in simplified model
Ke2Vals = KeP2Vals + 3600*pSol(9)./(1 + (pSol(10)*pSol(6)*3600./KeB2Vals).^pSol(8));

%% plot phase diagram with Hopf bifurcation (from running ddeBifCircadian.m)
fig3 = false;
supplFig = true;
YVals = logspace(-1,1.3,20);
MVals = logspace(-1,1.3,20);
[YMat, MMat] = meshgrid(YVals, MVals);
figure
if supplFig
    surf(YMat, MMat, circStoreCell{1}{2}*7.5, 'LineStyle','none','FaceColor','interp')
    cmap = loadInferno();
    cmap = cmap(81:end-80,:);
elseif fig3 %#ok<UNRCH>
    surf(YMat, MMat, circStoreCell{1}{1}/3600, 'LineStyle','none','FaceColor','interp')
    cmap = turbo(256);
    cmap = cmap(21:end-10,:);
end
daspect([1 1 1])
view([0 90])
hold on
plot3(BranchesStored{1}{1}, BranchesStored{1}{2}, 100*ones(size(BranchesStored{1}{1})),'k','LineWidth',3)
xlim([min(YMat(:)) max(YMat(:))])
ylim([min(MMat(:)) max(MMat(:))])
prettyGraph
colorbar
colormap(cmap)

xlabel('YAP/TAZ N/C')
ylabel('MRTF N/C')
set(gcf,'renderer','painters')
set(gca,'XScale','log')
set(gca,'YScale','log')

% paper figure with different treatment cases
colorOrder = linspecer(8);
colorOrder = vertcat([0, 0, 0], colorOrder);
if supplFig 
    stiffnessVec = [30, 0.3, 30, 30, 30, 30, 30, 30];
    CytDVec = [0, 0, 1, 0, 0, 0, 0, 0];
    LatBVec = [0, 0, 0, 0.2, 0, 0, 0, 0]; %technically LatA
    LATSVec = [1, 1, 1, 1, 3, 1, 1, 1]; % 1 is low density, 3 for high density
    blebbiVec = [0, 0, 0, 0, 0, 10, 0, 0];
    JasVec = [0,0,0,0,0,0,0,0];
    contactArea = [3000, 1000, 5000, 600, 1200, 4000, 1600, 900];
    markersVec = {'p','s','s','o','o','+','+','x'};
elseif fig3 %#ok<UNRCH>
    CytDVec = [0, 0, 0, 2, 5, 0, 0, 0];
    LatBVec = [0, 0, 0, 0, 0, 0, 0, 2];
    JasVec = [0, 0, 0, 0, 0, .1, .5, 0];
    stiffnessVec = [1e7, 10, 1, 1e7, 1e7, 1e7, 1e7, 1e7];
    LATSVec = [1, 1, 1, 1, 1, 1, 1, 1]; % 1 is low density, 3 for high density
    blebbiVec = [0, 0, 0, 0, 0, 0, 0, 0];
    contactArea = 3000*ones(size(CytDVec));
    markersVec = {'p','s','s','o','o','+','+','x'};
end
YAPTAZEq = zeros(size(stiffnessVec));
MRTFEq = zeros(size(stiffnessVec));
yStored = cell(size(stiffnessVec));
for i = 1:length(stiffnessVec)
    actinInhib =  1 / (1 + (LatBVec(i)/pSol(26))) + (1 + pSol(27))*JasVec(i) / pSol(28);
    cytoDConc = CytDVec(i);
    inhibVec = [actinInhib, 1, 1, 0, cytoDConc, LATSVec(i)]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)]
    inhibVec(7) = 1/(1 + blebbiVec(i)/1.0);
    inhibVec(8) = contactArea(i);
    inhibVec(9) = 1;
    [~,~, ~, ~, rawOutput] = conditionToOutputs(pSol,stiffnessVec(i),inhibVec,24*15*3600);
    t = rawOutput{1};
    y = rawOutput{2};
    ySS = rawOutput{3};
    YAPTAZEq(i) = ySS(15)/(ySS(17)+ySS(18));
    MRTFEq(i) = ySS(25)/ySS(26);
    plot3(YAPTAZEq(i),MRTFEq(i),100,markersVec{i},'MarkerSize',8,'LineWidth',2,'Color',colorOrder(i,:))
    yStored{i} = [t, y];
end
set(gcf,'renderer','painters')

%% plot dynamics for Fig 3
figure
groupsIdx = {1:3, [1,4,5], [1,6,7]};%, [1,8]};
refConc = 4.2;
legendStrings = {'Control', '10 kPa substrate', '1 kPa substrate',...
    '2 uM CytD', '5 uM CytD', '0.1 uM Jas', '0.5 uM Jas', '1 uM LatB'};
for i = 1:length(groupsIdx)
    subplot(length(groupsIdx),1,i)
    curGroupIdx = groupsIdx{i};
    for j = curGroupIdx
        oscDynamics = yStored{j}(:,3); %PER/CRY data
        t = yStored{j}(:,1)/(24*3600);
        [~,locs] = findpeaks(oscDynamics);
        tShift = t(locs(2));
        plot(t-tShift, refConc*yStored{j}(:,3), 'LineWidth', 1, 'Color', colorOrder(j,:)) % plot PER/CRY
        hold on
    end
    prettyGraph
    xlim([0 10])
    ylim([0 18])
    ylabel('Nuclear PER/CRY (nM)')
    xlabel('Time (days)')
    legend(legendStrings(curGroupIdx))
end
set(gcf,'renderer','painters')

%% plot all bifurcation diagrams
figure
spCount = 0;
if nCoupleTest 
    idx = 1:5;
elseif KdTest %#ok<UNRCH>
    idx = [2,4,6,7,8,9]
end
for i = idx
    spCount = spCount + 1;
    subplot(3,2,spCount)
    oscPeriod = circStoreCell{i}{1}/3600;
    oscDecayRate = circStoreCell{i}{3};
    oscPeriod(oscDecayRate<-.05) = nan;
    minPeriod = min(oscPeriod(:));
    oscPeriod(isnan(oscPeriod)) = .95*min(oscPeriod(:));
    surf(YMat, MMat, oscPeriod, 'LineStyle','none','FaceColor','interp')
    daspect([1 1 1])
    view([0 90])
    hold on
    plot3(BranchesStored{i}{1}, BranchesStored{i}{2}, 100*ones(size(BranchesStored{i}{1})),'k','LineWidth',2)
    xlim([min(YMat(:)) max(YMat(:))])
    ylim([min(MMat(:)) max(MMat(:))])
    prettyGraph
    colorbar
    % colormap spring
    cmap = turbo(256);
    % cmap = cmap + (1-cmap)*.15;
    cmap = cmap(21:end-10,:);
    cmap = vertcat([.7,.7,.7],cmap); %#ok<AGROW>
    colormap(cmap)
    
    xlabel('YAP/TAZ N/C')
    ylabel('MRTF N/C')
    set(gcf,'renderer','painters')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    if i == 6
        clim([.99*minPeriod, 31])
    else
        clim([.99*minPeriod, max(oscPeriod(:))])
    end
    clim([21.23 25.26])
end

%% test SS
pExp = 4;
tauB = p(1); pExpB = pExp; KeB = p(3);
KiB = p(4); KdBP = p(5); KdB = p(6);
tauP = p(7); pExpP0 = pExp; KeP = p(9);
KaP = p(10); KdP = p(11);
KeP1 = p(14); KiP = p(15);
KdR = p(16); tauR = p(17); 
pExpP1 = pExp;
KeB2 = 0.01;
KeR2 = 0.01;
KeP2 = 0:.001:.1;
PStar = zeros(size(KeP2));
for i = 1:length(KeP2)
    initGuess = [1, 1, 1];
    unstableSS = fsolve(@(ss) [KeB/(1+(ss(3)/KiB)^pExpB) + KeB2 - KdB*ss(1);
                               KeP/(1+(KaP/ss(1))^pExpP0) + KeP1*(1./(1+(ss(2)/KiP).^pExpP1)) + KeP2(i) - KdP*ss(2);
                               KeP/(1+(KaP/ss(1))^pExpP0) + KeP1*(1./(1+(ss(2)/KiP).^pExpP1)) + KeR2 - KdR*ss(3)],...
                initGuess);
    PStar(i) = unstableSS(2);
end
figure
plot(KeP2, PStar)
PStarApprox = (KeP1*KiP^pExp/KdP)^(1/(pExp+1)) + (1/((pExp+1)*KdP))*(KeP2 + KeP/(1+(KaP*KdB/KeB2)^pExp));
hold on
plot(KeP2, PStarApprox, '--')

% Functions defining the DDE system
function [t,y] = DDESolve2(timeSpan,p)
    tauB = p(1);    
    pExpB = p(2);
    KeB = p(3);
    KiB = p(4);
    % KdBP = p(5);
    KdB = p(6);
    tauP = p(7);
    pExpP0 = p(8);
    KeP = p(9);
    KaP = p(10);
    KdP = p(11);
    KeB2 = p(12);
    KeP2 = p(13);
    KeP1 = p(14);
    KiP = p(15);
    KdR = p(16);
    tauR = p(17);
    KeR2 = p(18);
    pExpP1 = p(19);
    
    fsolveOptions = optimoptions('fsolve','Display','off');
    ssVals = fsolve(@(ss) [KeB/(1+(ss(3)/KiB)^pExpB) + KeB2 - KdB*ss(1);
                  KeP/(1+(KaP/ss(1))^pExpP0) + KeP1/(1+(ss(2)/KiP)^pExpP1) + KeP2 - KdP*ss(2);
                  KeP/(1+(KaP/ss(1))^pExpP0) + KeP1/(1+(ss(2)/KiP)^pExpP1) + KeR2 - KdR*ss(3)],...
                [1,1,1], fsolveOptions);
    DDESol = dde23(@ddefun, [tauB,tauP,tauR], @history, timeSpan, ddeset('MaxStep',.1));
    t = DDESol.x';
    y = DDESol.y';

    function dy = ddefun(~,y,Z)
        RLag1 = Z(3,1);
        BLag2 = Z(1,2);
        PLag2 = Z(2,2);
        BLag3 = Z(1,3);
        PLag3 = Z(2,3);
        B = y(1);
        P = y(2);
        R = y(3);
        dy = [KeB/(1+(RLag1/KiB)^pExpB) + KeB2 - KdB*B;
              KeP/(1+(KaP/BLag2)^pExpP0) + KeP1/(1+(PLag2/KiP)^pExpP1) + KeP2 - KdP*P;
              KeP/(1+(KaP/BLag3)^pExpP0) + KeP1/(1+(PLag3/KiP)^pExpP1) + KeR2 - KdR*R];
    end
    
    function s = history(~)
        s = [5*ssVals(1); 0.2*ssVals(2); 1*ssVals(3)];
    end
end