%% two variables
% tauB1 = 8; pExpB = 2.5; KeB = 1;
% KiB = .04; KdBP = 0.5; KdB = 0.4;
% tauB2 = 8; pExpP = 2.5; KeP = 1;
% KaP = .5; KdP = .4;
% KeB2 = .1; KeP2 = .1;
% p = [tauB1, pExpB, KeB, KiB, KdBP, KdB, tauB2, pExpP, KeP, KaP, KdP, KeB2, KeP2];

p0 = [30112.1128026662/3600	2.49361433104013	0.00110195319176539*3600	0.0362906811897288	...
        0.000138718472336416*3600	4.01740850738952e-05*3600	28800/3600	1.13993216317370	0.000138008691550797*3600	...
        0.316725987150128	0.000225606126486809*3600	4.35158691215436e-05*3600    1.19993905179384e-05*3600]';
p = p0;
p(5) = 50;

timeSpan = [0 120];
KeB2Vals = 0:.05:1;
KeP2Vals = 0:.04:.8;
[KeB2Mat, KeP2Mat] = meshgrid(KeB2Vals, KeP2Vals);
% KeB2Mat = .1; KeP2Mat = .1;

period = zeros(size(KeB2Mat));
amplitude = zeros(size(KeB2Mat));
amplitudeRatio = zeros(size(KeB2Mat));
phaseLag = zeros(size(KeB2Mat));
timeFig = figure;

[numRows, numCols] = size(KeB2Mat);
parfor i = 1:numRows
    for j = 1:numCols
        pCur = p;
        figure(timeFig)
        pCur(12) = KeB2Mat(i,j);
        pCur(13) = KeP2Mat(i,j);
        [t2,y2] = DDESolve2(timeSpan,pCur);
        if any(isnan(y2(:,2)))
            y2(:,2) = 0;
        end
        if any(~isreal(y2(:,2)))
            y2(:,2) = real(y2(:,2));
        end
        if mod(i-1,5) == 0
%             subplot(3,1,j)
            plot(t2,y2(:,1),'LineWidth',1)
            prettyGraph
            xlabel('Time (hr)')
            ylabel('P')
            xlim([0 200])
            hold on
        end
        expressionFunc = pCur(9)./(1+(pCur(10)./y2(:,1)).^pCur(8));
        [pks,locs] = findpeaks(y2(:,2),'MinPeakProminence',.001);
        [troughs,troughLocs] = findpeaks(-y2(:,2),'MinPeakProminence',.001);
%         [pks2,locs2] = findpeaks(y2(:,1),'MinPeakProminence',.001);
%         [troughs2,troughLocs2] = findpeaks(-y2(:,1),'MinPeakProminence',.001);
        if length(pks)<=2 || length(troughs)<=2
            period(i,j) = nan;
            amplitude(i,j) = nan;
        else
            lastExtremum = min([length(pks),length(troughs)]); 
            amplitudeStored = pks(2:lastExtremum) - (-troughs(2:lastExtremum));
            if amplitudeStored(end) < 0.1*amplitudeStored(1) % cannot be decaying rapidly (not sustained osc)
                period(i,j) = nan;
                amplitude(i,j) = nan;
            else
                period(i,j) = mean(diff(t2(locs(2:end))));
                amplitude(i,j) = mean(amplitudeStored);
            end
        end
    end
end
figure
surf(KeB2Mat, KeP2Mat, period, 'LineStyle','none','FaceColor','interp')
daspect([.5 1 1])
view([0 90])

%% Parameter estimation for 2 DDE model
p0 = [8; 2.5; 1; .04; 0.5; 0.4; 8; 2.5; 1; .5; .4; 0; 0];
lowerLim = p0/2;
upperLim = p0*2;
% upperLim(end-1:end) = 1;
options = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon);
pSol = particleswarm(@pToObj_2Var, length(p0), lowerLim, upperLim,options);%p0);
pToObj_2Var(pSol)

%% 2 variables analytical (kinda)
% syms KeB KeP KdP KdBP KdB KiB KaP nB nP
nB = 2.5;
KeB = 1;
KiB = .04;
KdBP = .5;
KdBPVals = KdBP;% * logspace(-1,0,10);
KdB = 0.4;
nP = 2.5;
KeP = 1;
KaP = .5;
KdP = .4;
KeB2 = 0;% 0.1;
KeP2 = 0;%.1;
% KdPVals = .4*logspace(-1,1,100);
tauCur = 8;

testVals = KdB * logspace(-2,0,100);
tauCrit = zeros(size(KdBPVals));
periodEst = zeros(size(KdBPVals));
figure
hold on

tauB = tauCur;
tauP = tauCur;
p = {tauB, nB, KeB, KiB, KdBP, KdB, tauP, nP, KeP, KaP, KdP, KeB2, KeP2};


for k = 1:length(KdBPVals)
    KdBP = KdBPVals(k);
    p{5} = KdBP;
    for i = 1:length(testVals)
        KdB = testVals(i);
        p{13} = KeB2;
%         tauB = testVals(i);
%         p{1} = tauB;
%         % syms BStar PStar
%         solveOptions = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         SSEqns_no2 = @(star) [KeB/(1+(star(1)/KiB)^nB) - KdBP*star(1)*star(2) - KdB*star(1);...
%             KeP/(1+(KaP/star(1))^nP) - KdP*star(2)];
%         SSSol_ref = fsolve(SSEqns,[0,0],solveOptions);
%         SSEqns = @(star) [KeB/(1+(star(1)/KiB)^nB) - KdBP*star(1)*star(2) - KdB*star(1) + KeB2;...
%             KeP/(1+(KaP/star(1))^nP) - KdP*star(2) + KeP2];
%         SSSol = fsolve(SSEqns,[SSSol_ref(1),SSSol_ref(2)],solveOptions);
%         BStar = SSSol(1);
%         PStar = SSSol(2);
% 
%         phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
%         phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
% 
%         % phiModP = -KdBP*BStar*phiP/(omega^2 - KdP^2);
%         % define system of eqns
%         % eqn1 = sin(omega*tau) * (phiB - phiModP*KdP) - cos(omega*tau) * (phiModP*omega) == omega;
%         % eqn2 = sin(omega*tau) * (phiModP*omega) + cos(omega*tau) * (phiB - phiModP*KdP) == -KdBP*PStar - KdB;
%         % [omegaCrit,tauCrit] = solve([eqn1,eqn2], [omega,tau]);
% 
%         omegaTauEqns = @(omegaTau) [sin(omegaTau(1)*omegaTau(2)) * (phiB + KdP*KdBP*BStar*phiP/(omegaTau(1)^2 + KdP^2)) - cos(omegaTau(1)*omegaTau(2)) * ...
%             (-omegaTau(1)*KdBP*BStar*phiP/(omegaTau(1)^2 + KdP^2)) - omegaTau(1);...
%             sin(omegaTau(1)*omegaTau(2)) * (-omegaTau(1)*KdBP*BStar*phiP/(omegaTau(1)^2 + KdP^2)) + cos(omegaTau(1)*omegaTau(2)) * ...
%             (phiB + KdP*KdBP*BStar*phiP/(omegaTau(1)^2 + KdP^2)) + KdBP*PStar + KdB];
%         omegaTauSol = fsolve(omegaTauEqns,[.1,2]);
%         periodEst(i) = 2*pi*8 / (omegaTauSol(1)*omegaTauSol(2));
%         if ~isreal(omegaTauSol(2))
%             tauCrit(i) = abs(omegaTauSol(2));
%         else
%             tauCrit(i) = omegaTauSol(2);
%         end
        [x_lsq,x_lsqRed,omegaEst] = solve2Var(p);
        tauCrit_lsq(i) = x_lsq(4);
%         tauCrit_lsqRed(i) = x_lsqRed(4);
        omegaCrit_lsq(i) = x_lsq(3);
%         omegaCrit_lsqRed(i) = x_lsqRed(3);
        periodEst_lsq(i) = 2*pi*tauB / (x_lsq(3)*x_lsq(4));
%         periodEst_lsqRed(i) = 2*pi*tauB / (x_lsqRed(3)*x_lsqRed(4));
        amplRatio_lsq(i) = x_lsq(5);
        BStar(i) = x_lsq(1);
        PStar(i) = x_lsq(2);
%         amplRatio_lsqRed(i) = sqrt((x_lsqRed(3)^2 + KdP^2)/phiP^2);
%         periodShiftEst(i) = (1/(2*pi)) * (acos(KdP/sqrt(x_lsq(3)^2+KdP^2)) - x_lsq(4)*x_lsq(3));
%         resid1 = phiP*amplRatio_lsq(i)*cos(x_lsq(3)*x_lsq(4)) + KdP;
%         resid2 = phiP*amplRatio_lsq(i)*sin(x_lsq(3)*x_lsq(4)) - x_lsq(3);
%         residError(i) = resid1^2 + resid2^2;
%     end
%         BStar = x_lsqRed(1);
%         PStar = x_lsqRed(2);
%         phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
%         phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
%         omegaVals = 0:.01:10;
        omegaEq1 = @(omega) sin(omega.*tauB) .* (phiB + KdP*KdBP*BStar*phiP./(omega.^2 + KdP^2)) - cos(omega.*tauB) .* ...
            (-omega*KdBP*BStar*phiP./(omega.^2 + KdP^2)) - omega;
        omegaEq2 = @(omega) sin(omega*tauB) .* (-omega*KdBP*BStar*phiP./(omega.^2 + KdP^2)) + cos(omega*tauB) .* ...
            (phiB + KdP*KdBP*BStar*phiP./(omega.^2 + KdP^2)) + KdBP*PStar + KdB;
%         omegaEq3 = @(omega) cos(omega*tauB) - ...
%             (omega.^2 .* KdBP*BStar*phiP./(omega.^2+KdP^2) + (KdBP*PStar + KdB)*(phiB + KdP*KdBP*BStar*phiP./(omega.^2+KdP^2)))./...
%             ((phiB+KdP*(KdBP*BStar*phiP./(omega.^2+KdP^2))).^2 + (omega*KdBP*BStar*phiP./(omega.^2+KdP^2)).^2);
% %         options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');    
%         omegaEst = lsqnonlin(omegaEq1,omegaCrit_lsqRed(i)*.5,omegaCrit_lsqRed(i)*.1,omegaCrit_lsqRed(i));
%         periodEst_lsqRed(i) = 2*pi/omegaEst;
%         periodShiftEst(i) = (1/(2*pi)) * (acos(KdP/sqrt(omegaEst^2+KdP^2)) - tauB*omegaEst);
    end
    plot(testVals,omegaCrit_lsq)
    hold on
%     plot(testVals,omegaCrit_lsqRed)


end
% plot(KeP2Vals,periodEst)
% hold on
% plot(KeP2Vals,periodEst_lsq)

%% try new solving method for two DDEs
nB = 2.5;
KeB = 1;
KiB = .04;
KdBP = .5;
KdB = 0.4;
nP = 2.5;
KeP = 1;
KaP = .5;
KdP = .4;
KeB2 =  0.05;%0.1;
KeP2 = 0.05;%.1;
SSEqns = @(star) [KeB/(1+(star(1)/KiB)^nB) - KdBP*star(1)*star(2) - KdB*star(1) + KeB2;...
                  KeP/(1+(KaP/star(1))^nP) - KdBP*star(1)*star(2) - KdP*star(2) + KeP2];
SSSol = fsolve(SSEqns,[1; 1]);
BStar = SSSol(1);
PStar = SSSol(2);
phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
alphaLim = (phiB + phiP + KdP) / KdB;
% tauVec = 2:.1:10;
alphaVec = -0.8*alphaLim:.01:0.1*alphaLim;
omegaVec = zeros(size(alphaVec));
tauVec = zeros(size(alphaVec));
for i = 1:length(alphaVec)
    dynamicsEqns = @(dynVar) [dynVar(1)*(alphaVec(i)-1) + sin(dynVar(1)*dynVar(2))*(phiB+phiP);...
                              alphaVec(i)*KdP - KdB - cos(dynVar(1)*dynVar(2))*(phiB+phiP);
                              ];
    dynamicsSolve = fsolve(dynamicsEqns, [0.1; 8]);
    omegaVec(i) = dynamicsSolve(1);
%     alphaVec(i) = dynamicsSolve(2);
    tauVec(i) = dynamicsSolve(2);
end
figure
hold on
plot(alphaVec, omegaVec)

alphaPred1 = -(sin(omegaVec.*tauVec).*(phiP-phiB) - omegaVec) ./ omegaVec;
alphaPred2 = (cos(omegaVec.*tauVec).*(phiP-phiB) - KdB - 2*KdBP*PStar) ./ (KdP + 2*KdBP*BStar);
plot(alphaPred1, omegaVec)
plot(alphaPred2, omegaVec)
%% solve all
KeB2Vec = (0:.001:0.1);
KeP2Vec = 0:.001:0.1;
% KeB2Vec = .02*ones(size(KeP2Vec));
omegaVec = zeros(size(KeB2Vec));
tauVec = zeros(size(KeB2Vec));
alphaVec = zeros(size(KeB2Vec));
betaVec = zeros(size(KeB2Vec));
BStarVec = zeros(size(KeB2Vec));
PStarVec = zeros(size(KeB2Vec));
periodPred = zeros(size(KeB2Vec));
alphaPred = zeros(size(KeB2Vec));
betaPred = zeros(size(KeB2Vec));
for i = 1:length(KeB2Vec)
    SSEqns = @(star) [KeB/(1+(star(1)/KiB)^nB) - KdBP*star(1)*star(2) - KdB*star(1) + KeB2Vec(i);...
                      KeP/(1+(KaP/star(1))^nP) - KdBP*star(1)*star(2) - KdP*star(2) + KeP2Vec(i)];
    SSSol = fsolve(SSEqns,[1; 1]);
    BStar = SSSol(1);
    PStar = SSSol(2);
    phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
    phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
    dynamicsEqns = @(dynVar) [dynVar(1)*(dynVar(3)-1) + dynVar(4)*KdP + sin(dynVar(1)*dynVar(2))*(phiB+phiP);...
                              dynVar(3)*KdP - dynVar(4)*dynVar(1) - KdB - cos(dynVar(1)*dynVar(2))*(phiB+phiP);
                              dynVar(1) + dynVar(3) + sin(dynVar(1)*dynVar(2))*(phiP-phiB) + KdP*dynVar(4) + 2*KdBP*dynVar(4)*BStar;
                              -dynVar(4)*dynVar(1) - cos(dynVar(1)*dynVar(2))*(phiP-phiB) - KdB - KdP*dynVar(3) - 2*KdBP*(BStar*dynVar(3)+PStar)];
%     options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%     dynamicsSolve = fsolve(dynamicsEqns, [0.1; 8; -1; 0]);%,options);
    dynamicsSolve = lsqnonlin(dynamicsEqns, [0.1; 8; -1; 0]);
    omegaVec(i) = dynamicsSolve(1);
    tauVec(i) = dynamicsSolve(2);
    alphaVec(i) = dynamicsSolve(3);
    betaVec(i) = dynamicsSolve(4);
    BStarVec(i) = BStar;
    PStarVec(i) = PStar;
    dynamicsSolve2 = lsqnonlin(dynamicsEqns, [0.1; 8; -1; 0], [0; 8; -10; -pi], [1; 8; -.1; pi]);
    periodPred(i) = 2*pi/dynamicsSolve2(1);
    alphaPred(i) = dynamicsSolve2(3);
    betaPred(i) = dynamicsSolve2(4);
end
figure
plot(KeP2Vec,periodPred)

%% solve for omega and alpha near hopf
p = [pSol(1:12),pSol(15)];
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
KdBP = 0;
KdBVals = 0:.01:1;
KeBVals = KeB*[.1,.2,.5,1,2,5,10];
bifFig = figure;
hold on
omegaFig = figure;
hold on
for j = 1:length(KeBVals)
    KeB = KeBVals(j);
    KeB2Pred = zeros(size(KdBVals));
    omegaPred = zeros(size(KdBVals));
    for i = 1:length(KdBVals)
        KdB = KdBVals(i);
        SSEqns = @(star) [KeB/(1+(star(1)/KiB)^nB) - KdBP*star(1)*star(2) - KdB*star(1) + KeB2;...
            KeP/(1+(KaP/star(1))^nP) - KdBP*star(1)*star(2) - KdP*star(2) + KeP2];
        simplEqns = @(var) [KeB/(1+(var(1)/KiB)^nB) - KdB*var(1) + var(3);
            (KeB*(nB/KiB)*(var(1)/KiB)^(nB-1) / (1 + (var(1)/KiB)^nB)^2)*cos(var(2)*tau) + KdB;
            (KeB*(nB/KiB)*(var(1)/KiB)^(nB-1) / (1 + (var(1)/KiB)^nB)^2)*sin(var(2)*tau) - var(2)];
        simplSol = lsqnonlin(simplEqns,[0.1; 0.25; 0.06], [0; 0; 0], [1; 1; 1]);
        % BOnly = fsolve(@(star) KeB/(1+(star/KiB)^nB) - KdB*star + KeB2, 1);
        % BStar = SSSol(1);
        % PStar = SSSol(2);
        BStar = simplSol(1);
        phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
        phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
        % dynEqns = @(dynVar) [(KdP*dynVar(2)-KdB)^2 + (dynVar(1)*(dynVar(2)-1))^2 - (phiP + phiB)^2;...
        %     (2*KdBP*(BStar*dynVar(2)+PStar)-KdP*dynVar(2)-KdB)^2 + (dynVar(1)*(dynVar(2)+1))^2 - (phiP-phiB)^2];
        % dynSolve = lsqnonlin(dynEqns,[.1;-1],[0;-20],[1;-1]);
        KeB2Pred(i) = simplSol(3);
        omegaPred(i) = sqrt(phiB^2 - KdB^2);
        PEqns = @(pPar) [-omegaPred*pPar(1)*sin(pPar(2)) - cos(omegaPred*tau)*phiP + KdP*pPar(1)*cos(pPar(2));
            omegaPred*pPar(1)*cos(pPar(2)) + sin(omegaPred*tau)*phiP + KdP*pPar(1)*sin(pPar(2))];
        PSol = lsqnonlin(PEqns, [2.1966; pi], [0; 0], [10; 2*pi]);
    end
    figure(bifFig)
    plot(KdBVals,KeB2Pred)
    figure(omegaFig)
    plot(KeB2Pred, 2*pi./omegaPred)
end

%% treat KeB2 as unknown
% tauB = 8;
% tauP = 8;
p = [pSol(1:12),pSol(15)];
% p = {tauB, nB, KeB, KiB, KdBP, KdB, tauP, nP, KeP, KaP, KdP, KeB2, KeP2};
KdBPVec = 0:.01:2;%0:.001:1;
p(13) = p(13)/5;
KeB2Pred = zeros(size(KdBPVec));
% KeP2Pred = zeros(size(KdBPVec));
omegaPred = zeros(size(KdBPVec));
alphaPred = zeros(size(KdBPVec));
var0 = [0.1029, 2, .0610, .2611, 2.1966, 2.1862]';
% var0 = [.0914, 1.6768, .0877, .2628, 3.8046, 2.1474];
for i = 1:length(KdBPVec)
    p(5) = KdBPVec(i)/3600;
    curPred = hopfAnalysis(p, var0);
    KeB2Pred(i) = curPred(3);
%     KeP2Pred(i) = curPred(4);
    omegaPred(i) = curPred(4);
    alphaPred(i) = curPred(5);
    var0 = curPred;
end
phiB = p(3)*3600*(p(2)/p(4))*(curPred(1)/p(4))^(p(2)-1) / (1 + (curPred(1)/p(4))^p(2))^2;
phiP = p(9)*3600*p(8)*(p(10)^p(8)/curPred(1)^(p(8)+1)) / (1 + (p(10)/curPred(1))^p(8))^2;


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
    fullEqns(varPred)

    function F = fullEqns(var)
        % var = [BStar,PStar,KeB2,KeP2,omega,alpha]
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


function [x,xRed,omegaEst] = solve2Var(p)
    nB = p{2};
    KeB = p{3};
    KiB = p{4};
    KdBP = p{5};
    KdB = p{6};
    nP = p{8};
    KeP = p{9};
    KaP = p{10};
    KdP = p{11};
    KeB2 = p{12};
    KeP2 = p{13};
    tauB = p{1};
    xSS = fsolve(@findSS,[.1,.1]);
    BStar = xSS(1);
    PStar = xSS(2);
    options = optimoptions('particleswarm','HybridFcn','fmincon');%'Algorithm','trust-region-reflective');%,'TolFun',1e-3);
%     x = lsqnonlin(@leastSq2Var,[.2,2,1,pi],[.1,1,0,0],[1,100,100,2*pi],options);
%     x = particleswarm(@leastSq2Var,4,[0,0,0,0],[10,100,100,2*pi],options);
%     x = lsqnonlin(@leastSq2VarRed,[.1,2],[0,1],[100,100],options);
%     x = [x,-1];
    x = fsolve(@leastSq2Var,[.1,2,1,pi]);%,options);
    x = [xSS,x];
    xRed = [];
    omegaEst = [];
%     xRed = nan;
%     xRed = lsqnonlin(@leastSq2VarRed,[0.1,0.1,.1,2],zeros(1,4),[100*ones(1,4)],options);
%     BStar = xRed(1);
%     PStar = xRed(2);
%     omegaEst = lsqnonlin(@leastSqOmega,.5*xRed(3),0.2*xRed(3),0.8*xRed(3));

    function F = leastSq2Var(x)

%         BStar = x(1);
%         PStar = x(2);
        omega = x(1);
        tau = x(2);
        alpha = x(3);
        beta = x(4);
%         KdB = x(6);
%         KdP = x(7);

        phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
        phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
% 
%         F(1) = KeB/(1+(BStar/KiB)^nB) - KdBP*BStar*PStar - KdB*BStar + KeB2;
%         F(2) = KeP/(1+(KaP/BStar)^nP) - KdBP*BStar*PStar - KdP*PStar + KeP2;
        F(1) = -omega + phiB*sin(omega*tau)-KdBP*(BStar*alpha*sin(beta));
        F(2) = -phiB*cos(omega*tau) - KdBP*(BStar*alpha*cos(beta)+PStar)-KdB;
        F(3) = -omega - phiP*sin(omega*tau)+KdBP*(PStar*sin(beta)/alpha);
        F(4) = phiP*cos(omega*tau) - KdBP*(BStar + PStar*cos(beta)/alpha) - KdP;
    end

%     function F = leastSq2VarRed(x)
% 
%         BStar = x(1);
%         PStar = x(2);
%         omega = x(3);
%         tau = x(4);
%         %         amplRatio = x(5);
%         %         KdB = x(6);
%         %         KdP = x(7);
% 
%         phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
%         phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
% 
%         F(1) = KeB/(1+(BStar/KiB)^nB) - KdBP*BStar*PStar - KdB*BStar + KeB2;
%         F(2) = KeP/(1+(KaP/BStar)^nP) - KdP*PStar + KeP2;
%         F(3) = sin(omega*tau) * (phiB + KdP*KdBP*BStar*phiP/(omega^2 + KdP^2)) - cos(omega*tau) * ...
%             (-omega*KdBP*BStar*phiP/(omega^2 + KdP^2)) - omega;
%         F(4) = sin(omega*tau) * (-omega*KdBP*BStar*phiP/(omega^2 + KdP^2)) + cos(omega*tau) * ...
%             (phiB + KdP*KdBP*BStar*phiP/(omega^2 + KdP^2)) + KdBP*PStar + KdB;
%     end
% 
%     function F = leastSqOmega(omega)
% 
%         phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
%         phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;
% 
%         F(1) = sin(omega*tauB) * (phiB + KdP*KdBP*BStar*phiP/(omega^2 + KdP^2)) - cos(omega*tauB) * ...
%             (-omega*KdBP*BStar*phiP/(omega^2 + KdP^2)) - omega;
%         F(2) = sin(omega*tauB) * (-omega*KdBP*BStar*phiP/(omega^2 + KdP^2)) + cos(omega*tauB) * ...
%             (phiB + KdP*KdBP*BStar*phiP/(omega^2 + KdP^2)) + KdBP*PStar + KdB;
%     end

    function F = findSS(x)
        BStar = x(1);
        PStar = x(2);

        phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
        phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;

        F(1) = KeB/(1+(BStar/KiB)^nB) - KdBP*BStar*PStar - KdB*BStar + KeB2;
        F(2) = KeP/(1+(KaP/BStar)^nP) - KdBP*BStar*PStar - KdP*PStar + KeP2;
    end

    function F = leastSq2VarRed(x)
        omega = x(1);
        tau = x(2);
        A = -1;

        phiB = KeB*(nB/KiB)*(BStar/KiB)^(nB-1) / (1 + (BStar/KiB)^nB)^2;
        phiP = KeP*nP*(KaP^nP/BStar^(nP+1)) / (1 + (KaP/BStar)^nP)^2;

        F(1) = sin(omega*tau) * (phiB+A*phiP) - omega*(1-A);
        F(2) = -cos(omega*tau) * (phiB+A*phiP) - KdB + A*KdP;
    end
end

function [t,y] = DDESolve2(timeSpan,p)
    tauB1 = p(1);    
    pExpB = p(2);
    KeB = p(3);
    KiB = p(4);
    KdBP = p(5);
    KdB = p(6);
    tauB2 = p(7);
    pExpP = p(8);
    KeP = p(9);
    KaP = p(10);
    KdP = p(11);
    KeB2 = p(12);
    KeP2 = p(13);

    ssVals = fsolve(@(ss) [KeB/(1+(ss(1)/KiB)^pExpB) + KeB2 - KdBP*ss(1)*ss(2) - KdB*ss(1);
                           KeP/(1+(KaP/ss(1))^pExpP) + KeP2 - KdBP*ss(1)*ss(2) - KdP*ss(2)],...
                [0.1,0.1]);
    DDESol = dde23(@ddefun, [tauB1,tauB2], @history, timeSpan, ddeset('MaxStep',.1));
    t = DDESol.x';
    y = DDESol.y';

    function dy = ddefun(t,y,Z)
        BLag1 = Z(1,1);
        BLag2 = Z(1,2);
        B = y(1);
        P = y(2);
%         KeB2Cur = KeB2 * (1 - t/timeSpan(2));
        dy = [KeB/(1+(BLag1/KiB)^pExpB) + KeB2 - KdBP*B*P - KdB*B;
              KeP/(1+(KaP/BLag1)^pExpP) + KeP2 - KdBP*B*P - KdP*P];
    end
    
    function s = history(t)
        s = [5*ssVals(1); 0.2*ssVals(2)];
    end
end

function obj = pToObj_2Var(p)
    tauB1 = p(1);    
    pExpB = p(2);
    KeB = p(3);
    KiB = p(4);
    KdBP = p(5);
    KdB = p(6);
    tauB2 = p(7);
    pExpP = p(8);
    KeP = p(9);
    KaP = p(10);
    KdP = p(11);
    KeB2 = p(12);
    KeP2 = p(13);
    timeSpan = [0 500];

    DDESol = dde23(@ddefun, [tauB1,tauB2], @history, timeSpan, ddeset('MaxStep',.1));
    t = DDESol.x';
    y = DDESol.y';
    [pks,locs] = findpeaks(y(:,2),'MinPeakProminence',.1);
    [troughs,~] = findpeaks(-y(:,2),'MinPeakProminence',.1);
    if length(pks)>10 && length(troughs)>10
        period = mean(diff(t(locs(2:end))));
        lastExtremum = min([length(pks),length(troughs)]); 
        amplitudeVec = pks(2:lastExtremum) - (-troughs(2:lastExtremum));
        if amplitudeVec(end) < 0.8*amplitudeVec(1) % cannot be decaying rapidly (not sustained osc)
            period = t(end);
        end
%         amplitude = mean(pks(2:end)) - mean(-troughs(2:end));
    else
        period = t(end);
%         amplitude = mean(y(:,2));
    end
    obj = (period - 24)^2;


    function dy = ddefun(t,y,Z)
        BLag1 = Z(1,1);
        BLag2 = Z(1,2);
        B = y(1);
        P = y(2);
        dy = [KeB/(1+(BLag1/KiB)^pExpB) + KeB2 - KdBP*B*P - KdB*B;
              KeP/(1+(KaP/BLag2)^pExpP) + KeP2 - KdBP*B*P - KdP*P];
    end
    
    function s = history(t)
        s = zeros(2,1);
    end
end