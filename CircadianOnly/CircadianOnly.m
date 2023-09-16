%% two variables
% tauB1 = 8; pExpB = 2.5; KeB = 1;
% KiB = .04; KdBP = 0.5; KdB = 0.4;
% tauB2 = 8; pExpP = 2.5; KeP = 1;
% KaP = .5; KdP = .4;
% KeB2 = .1; KeP2 = .1;
% p = [tauB1, pExpB, KeB, KiB, KdBP, KdB, tauB2, pExpP, KeP, KaP, KdP, KeB2, KeP2];

% p0 = [30112.1128026662/3600	2.49361433104013	0.00110195319176539*3600	0.0362906811897288	...
%         0.000138718472336416*3600	4.01740850738952e-05*3600	28800/3600	1.13993216317370	0.000138008691550797*3600	...
%         0.316725987150128	0.000225606126486809*3600	4.35158691215436e-05*3600    1.19993905179384e-05*3600]';
% p = p0;
% p(5) = 0.5;
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 1/3600; .5; 0.4/3600;...
0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
100; 1; 10; 2; 2];
p = p0(1:13);%pSol(1:13);

% convert to units of hours
% p(12:13) = 0;
p([1,7]) = p([1,7])/3600;
p([3,5,6,9,11,12,13]) = p([3,5,6,9,11,12,13])*3600;

timeSpan = [0 960];
YVals = logspace(-1,1.3,40);
MVals = logspace(-1,1.3,40);
[YMat, MMat] = meshgrid(YVals, MVals);
if ~exist('pSol','var')
    error('Load in mechano-Circadian pSol vector first')
end
CytoConv = 1.3851e6;
NucConv = 3.3122e5;
Ytot = 1.4784e6;
Mtot = 1e6;
YConcMat = YMat*Ytot./(CytoConv + YMat*NucConv);
MConcMat = MMat*Mtot./(CytoConv + MMat*NucConv);
KeB2Mat = 3600*(pSol(12)*YConcMat.^2./(pSol(13)^2+YConcMat.^2) +...
                pSol(23)*MConcMat.^2./(pSol(24)^2+MConcMat.^2));
KeP2Mat = 3600*(pSol(15)*MConcMat.^2./(pSol(16)^2+MConcMat.^2) +...
                pSol(20)*YConcMat.^2./(pSol(21)^2+YConcMat.^2)); 
% KeB2Vals = .03;%0:.005:.1;
% KeP2Vals = 0;%0:.04:.8;
% [KeB2Mat, KeP2Mat] = meshgrid(KeB2Vals, KeP2Vals);
% % KeB2Mat = .1; KeP2Mat = .1;

KdBPVals = [1,0,.1,.3,1,1,1,1,1,1,...
            1,1,1,1,1,1,1]*p(5); % KdBP (p(5))
KdBVals = [1,1,1,1,.1,.3,3,1,1,1,...
           1,1,1,1,1,1,1]*p(6); % KdB (p(6))
KdPVals = [1,1,1,1,1,1,1,3,10,30,...
           1,1,1,1,1,1,1]*p(11); % KdP (p(11))
tauBVals = [1,1,1,1,1,1,1,1,1,1,...
            .1,.3,3,1,1,1,1]*p(1); % tauB (p(1))
tauPVals = [1,1,1,1,1,1,1,1,1,1,...
            1,1,1,.1,.3,3,10]*p(7); % tauP (p(7))
circStoreCell = cell(size(KdBPVals));

for k = 2:10%length(KdBPVals)
    p(5) = KdBPVals(k); p(6) = KdBVals(k); p(11) = KdPVals(k);
    p(1) = tauBVals(k); p(7) = tauPVals(k);
    period = zeros(size(KeB2Mat));
    oscDecayRate = zeros(size(KeB2Mat));
    amplitude = zeros(size(KeB2Mat));
    amplitudeRatio = zeros(size(KeB2Mat));
    phaseLag = zeros(size(KeB2Mat));
    % timeFig = figure;

    [numRows, numCols] = size(KeB2Mat);
    for i = 1:numRows
        fprintf('Test %d, row %d of %d\n', k, i, numRows)
        for j = 1:numCols
            pCur = p;
            % figure(timeFig)
            pCur(12) = KeB2Mat(i,j);
            pCur(13) = KeP2Mat(i,j);
            [t2,y2] = DDESolve2(timeSpan,pCur);
            if any(isnan(y2(:,2)))
                y2(:,2) = 0;
            end
            if any(~isreal(y2(:,2)))
                y2(:,2) = real(y2(:,2));
            end
            %         if mod(i-1,5) == 0
            % %             subplot(3,1,j)
            %             plot(t2,y2(:,1),'LineWidth',1)
            %             prettyGraph
            %             xlabel('Time (hr)')
            %             ylabel('P')
            %             xlim([0 200])
            %             hold on
            %         end
            expressionFunc = pCur(9)./(1+(pCur(10)./y2(:,1)).^pCur(8));
            [period(i,j), amplitude(i,j), oscDecayRate(i,j)] = circOscAnalysis(t2*3600, y2(:,2));
        end
    end
    circStoreCell{k} = {period, amplitude, oscDecayRate};
end

%%
figure
for i = 1%:10%length(circStoreCell)
    % subplot(3,3,i)
    oscPeriod = circStoreCell{i}{1}/3600;
    oscDecayRate = circStoreCell{i}{3};
    oscPeriod(oscDecayRate<-.01) = nan;
    % for j = 1:size(oscPeriod,1)
    %     curRow = oscPeriod(j,:);
    %     lastPeriodIdx = find(~isnan(curRow),1,'last');
    %     curRow(isnan(curRow)) = 0;%curRow(lastPeriodIdx);
    %     oscPeriod(j,:) = curRow;
    % end
    surf(YMat, MMat, oscPeriod, 'LineStyle','none','FaceColor','interp')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    daspect([1 1 1])
    view([0 90])
    drawnow
    colorbar
    % clim([-.01 0])
    % clim([32 33])
    xlim([.1 20])
    ylim([.1 20])
    clim([min(oscPeriod(oscPeriod>0)) max(oscPeriod(:))])
end

%% plot with Hopf
YVals = logspace(-1,1.3,40);
MVals = logspace(-1,1.3,40);
[YMat, MMat] = meshgrid(YVals, MVals);
figure
surf(YMat, MMat, circStoreCell{1}{1}/3600, 'LineStyle','none','FaceColor','interp')
daspect([1 1 1])
view([0 90])
hold on
plot3(BranchesStored{1}{1}, BranchesStored{1}{2}, 100*ones(size(BranchesStored{1}{1})),'k','LineWidth',3)
xlim([min(YMat(:)) max(YMat(:))])
ylim([min(MMat(:)) max(MMat(:))])
% clim([24 28])
prettyGraph
colorbar
% colormap spring
cmap = turbo(256);
cmap = cmap(21:end-10,:);
% cmap = cmap + (1-cmap)*.15;
colormap(cmap)

xlabel('YAP/TAZ N/C')
ylabel('MRTF N/C')
set(gcf,'renderer','painters')
set(gca,'XScale','log')
set(gca,'YScale','log')

% paper figure with different treatment cases
% turboMap = parula;
colorOrder = linspecer(6);
colorOrder = vertcat([0, 0, 0], colorOrder);
CytDVec = [0, 0, 0, 2, 5, 0, 0];
JasVec = [0, 0, 0, 0, 0, .1, .5];
stiffnessVec = [1e7, 10, 1, 1e7, 1e7, 1e7, 1e7];
markersVec = {'p','s','s','o','o','+','+'};
YAPTAZEq = zeros(size(stiffnessVec));
MRTFEq = zeros(size(stiffnessVec));
yStored = cell(size(stiffnessVec));
for i = 1:length(stiffnessVec)
    actinInhib = 1 + (1 + pSol(27))*JasVec(i) / pSol(28);
    cytoDConc = CytDVec(i);
    inhibVec = [actinInhib, 1, 1, 1, cytoDConc, 1]; %[actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)] 
    [~,~, ~, ~, rawOutput] = conditionToOutputs(pSol,stiffnessVec(i),inhibVec,24*15*3600);
    t = rawOutput{1};
    y = rawOutput{2};
    ySS = rawOutput{3};
    YAPTAZEq(i) = ySS(15)/(ySS(17)+ySS(18));
    MRTFEq(i) = ySS(25)/ySS(26);
    plot3(YAPTAZEq(i),MRTFEq(i),100,markersVec{i},'MarkerSize',8,'LineWidth',2,'Color',colorOrder(i,:))
    yStored{i} = [t, y];
end
%%
figure
groupsIdx = {1:3, [1,4,5], [1,6,7]};
refConc = 13.5;
legendStrings = {'Control', '10 kPa substrate', '1 kPa substrate',...
    '2 uM CytD', '5 uM CytD', '0.1 uM Jas', '0.5 uM Jas'};
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
    ylim([0 22])
    ylabel('Nuclear PER/CRY (nM)')
    xlabel('Time (days)')
    legend(legendStrings(curGroupIdx))
end
set(gcf,'renderer','painters')

%% control bifurcation, both period and ampl
YVals = logspace(-1,1.3,40);
MVals = logspace(-1,1.3,40);
[YMat, MMat] = meshgrid(YVals, MVals);
figure
conversionVec = [1/3600, 1, 1];
for i = 1:3
    subplot(3,1,i)
    surf(YMat, MMat, circStoreCell{1}{i}*conversionVec(i), 'LineStyle','none','FaceColor','interp')
    daspect([1 1 1])
    view([0 90])
    hold on
    plot3(BranchesStored{1}{1}, BranchesStored{1}{2}, 100*ones(size(BranchesStored{1}{1})),'k','LineWidth',3)
    xlim([min(YMat(:)) max(YMat(:))])
    ylim([min(MMat(:)) max(MMat(:))])
    % clim([24 28])
    prettyGraph
    colorbar
    cmap = turbo(256);
    cmap = cmap(21:end-10,:);
    colormap(cmap)
    xlabel('YAP/TAZ N/C')
    ylabel('MRTF N/C')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
end
set(gcf,'renderer','painters')

%% plot all bifurcation diagrams
figure
spCount = 0;
for i = [2,4,5,7,8,10]
    spCount = spCount + 1;
    subplot(3,2,spCount)
    oscPeriod = circStoreCell{i}{1}/3600;
    oscDecayRate = circStoreCell{i}{3};
    oscPeriod(oscDecayRate<-.01) = nan;
    minPeriod = min(oscPeriod(:));
    oscPeriod(isnan(oscPeriod)) = .95*min(oscPeriod(:));
    surf(YMat, MMat, oscPeriod, 'LineStyle','none','FaceColor','interp')
    daspect([1 1 1])
    view([0 90])
    hold on
    plot3(BranchesStored{i}{1}, BranchesStored{i}{2}, 100*ones(size(BranchesStored{i}{1})),'k','LineWidth',2)
    xlim([min(YMat(:)) max(YMat(:))])
    ylim([min(MMat(:)) max(MMat(:))])
    % clim([25 35])
    prettyGraph
    colorbar
    % colormap spring
    cmap = turbo(256);
    % cmap = cmap + (1-cmap)*.15;
    cmap = cmap(21:end-10,:);
    cmap = vertcat([.7,.7,.7],cmap);
    colormap(cmap)
    
    xlabel('YAP/TAZ N/C')
    ylabel('MRTF N/C')
    set(gcf,'renderer','painters')
    set(gca,'XScale','log')
    set(gca,'YScale','log')

    clim([.99*minPeriod, max(oscPeriod(:))])
end


%% solve for Hopf bifurcation, treating KeB2 or KeP2 as unknown
p = [pSol(1:12),pSol(15)];
% p = {tauB, nB, KeB, KiB, KdBP, KdB, tauP, nP, KeP, KaP, KdP, KeB2, KeP2};
KeP2Vec = 0:.01:0.8;%0:.001:1;
p(5) = 0.5/3600;
KeB2Pred = zeros(size(KeP2Vec));
% KeP2Pred = zeros(size(KdBPVec));
omegaPred = zeros(size(KeP2Vec));
alphaPred = zeros(size(KeP2Vec));
var0 = [0.2456, 0.2274, .04, .2611, 0.5, pi]';
% var0 = [.0914, 1.6768, .0877, .2628, 3.8046, 2.1474];
for i = 1:length(KeP2Vec)
    p(13) = KeP2Vec(i)/3600;
    curPred = hopfAnalysis(p, var0);
    KeB2Pred(i) = curPred(3);
    omegaPred(i) = curPred(4);
    alphaPred(i) = curPred(5);
    var0 = curPred;
end
phiB = p(3)*3600*(p(2)/p(4))*(curPred(1)/p(4))^(p(2)-1) / (1 + (curPred(1)/p(4))^p(2))^2;
phiP = p(9)*3600*p(8)*(p(10)^p(8)/curPred(1)^(p(8)+1)) / (1 + (p(10)/curPred(1))^p(8))^2;
figure
plot(KeB2Pred, KeP2Vec)

%% solve simplified case of KdBP=0
KeB = 1;
KiB = 0.04;
KdB = 0.4;
KeB2Vec = 0:.001:.1;
BStarVec = zeros(size(KeB2Vec));
phiBVec = zeros(size(KeB2Vec));
for i = 1:length(KeB2Vec)
    KeB2 = KeB2Vec(i);
    curFun = @(BStar) KeB/(1+(BStar/KiB)^2) + KeB2 - KdB*BStar;
    BStarVec(i) = fsolve(curFun, 0.1);
    BHat = BStarVec(i)/KiB;
    phiBVec(i) = (2*KeB/KiB) * (BHat/(1+BHat^2)^2);
end


function varPred = hopfAnalysis(p, var0)
    nB = p(2);
    KeB = p(3)*3600;
    KiB = p(4);
    KdBP = p(5)*3600;
    KdB = p(6)*3600;
    nP = p(8);
    KeP = p(9)*3600;
    KaP = p(10);
    KdP = p(11)*3600;
    KeB2 = p(12)*3600;
    KeP2 = p(13)*3600;
    tau = p(1)/3600;
    varPred = lsqnonlin(@fullEqns, var0, [0, 0, 0, 0, 0, 0]', [10, 10 , 2, 1, 10, 2*pi]');
%     fullEqns(varPred)

    function F = fullEqns(var)
        % var = [BStar,PStar,KeB2,omega,alpha,beta]
        BStar = var(1);
        PStar = var(2);
        KeB2 = var(3);
%         KeP2 = var(4);
        omega = var(4);
        alpha = var(5);
        beta = var(6);
        phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
        phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
%         FSimpl = [KeB/(1+(BStar/KiB)^nB) - KdBP*BStar*PStar - KdB*BStar + KeB2;
%                   KeP/(1+(KaP/BStar)^nP) - KdBP*BStar*PStar - KdP*PStar + KeP2;                     
%                   phiB*cos(omega*tau) + KdB;
%                   phiB*sin(omega*tau) - omega;
%                   -omega*alpha*sin(beta) - cos(omega*tau)*phiP + KdP*alpha*cos(beta);
%                   omega*alpha*cos(beta) + sin(omega*tau)*phiP + KdP*alpha*sin(beta)];
        F = [KeB/(1+(BStar/KiB)^nB) - KdBP*BStar*PStar - KdB*BStar + KeB2;
             KeP/(1+(KaP/BStar)^nP) - KdBP*BStar*PStar - KdP*PStar + KeP2;
             phiB*cos(omega*tau) + KdBP*(BStar*alpha*cos(beta) + PStar) + KdB;
             omega - phiB*sin(omega*tau) + KdBP*BStar*alpha*sin(beta);
             -(phiP/alpha)*cos(beta+omega*tau) + KdBP*(BStar + PStar*cos(beta)/alpha) + KdP;
             omega + (phiP/alpha)*sin(beta+omega*tau) - KdBP*PStar*sin(beta)/alpha];
    end
end


function [t,y] = DDESolve2(timeSpan,p)
    tauB = p(1);    
    pExpB = p(2);
    KeB = p(3);
    KiB = p(4);
    KdBP = p(5);
    KdB = p(6);
    tauP = p(7);
    pExpP = p(8);
    KeP = p(9);
    KaP = p(10);
    KdP = p(11);
    KeB2 = p(12);
    KeP2 = p(13);
    
    fsolveOptions = optimoptions('fsolve','Display','off');
    ssVals = fsolve(@(ss) [KeB/(1+(ss(1)/KiB)^pExpB) + KeB2 - KdBP*ss(1)*ss(2) - KdB*ss(1);
                           KeP/(1+(KaP/ss(1))^pExpP) + KeP2 - KdBP*ss(1)*ss(2) - KdP*ss(2)],...
                [0.1,0.1], fsolveOptions);
    DDESol = dde23(@ddefun, [tauB,tauP], @history, timeSpan, ddeset('MaxStep',.1));
    t = DDESol.x';
    y = DDESol.y';

    function dy = ddefun(t,y,Z)
        BLag = Z(1,1);
        PLag = Z(1,2);
        B = y(1);
        P = y(2);
%         KeB2Cur = KeB2 * (1 - t/timeSpan(2));
        dy = [KeB/(1+(BLag/KiB)^pExpB) + KeB2 - KdBP*B*P - KdB*B;
              KeP/(1+(KaP/PLag)^pExpP) + KeP2 - KdBP*B*P - KdP*P];
    end
    
    function s = history(t)
        s = [5*ssVals(1); 0.2*ssVals(2)];
    end
end