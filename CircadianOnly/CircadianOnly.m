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
p(5) = 0.5;

timeSpan = [0 240];
KeB2Vals = .03;%0:.005:.1;
KeP2Vals = 0;%0:.04:.8;
[KeB2Mat, KeP2Mat] = meshgrid(KeB2Vals, KeP2Vals);
% KeB2Mat = .1; KeP2Mat = .1;

period = zeros(size(KeB2Mat));
amplitude = zeros(size(KeB2Mat));
amplitudeRatio = zeros(size(KeB2Mat));
phaseLag = zeros(size(KeB2Mat));
timeFig = figure;

[numRows, numCols] = size(KeB2Mat);
for i = 1:numRows
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