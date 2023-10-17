function [t,y, SSVals] = MechanoCircadianModel(timeSpan, stiffnessParam, pSol, inhibVec, varargin)
% Function implementing DDE system for modeling
% mechanotransduction-Circadian coupling. Calls MechanoSS to find SS
% nuclear concentrations of YAP/TAZ and MRTF.
% Emmet Francis, 2023
%
% input:
%     * timeSpan: two element vector [tStart tEnd] in seconds
%
%     * stiffnessParam is a 2 element vector with the stiffness in kPa (first
%       element) and the time scale of stiffness increase in s (second el)
%       only stiffnessParam(1) is used here (stiffness is assumed to be
%       equal to this value at SS)
%
%     * pSol: parameter solution vector - length 34, contains all calibrated
%       values from mechano-Circadian model (see Tables 1-2 in paper). Note
%       that nB, nP, nYB, nMP, nYP, and nMB are all Hill coefficients that were
%       not included in parameter estimation, but were fixed at 2.0.
%       In order, the parameters are:
%       {'tauB','nB','KeB0','KiB','KdBP','KdB',...
%       'tauP','nP','KeP0','KaP','KdP',...
%       'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
%       'ROCKInhibSens', 'StiffThresh',...
%       'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
%       'LatBSens','JasMag','JasSens','KdLuc',...
%       'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim','Kcap'};
%
%     * inhibVec: vector of inhibition parameters (length=9)
%           1: actin polym inhibition: factor multiplying kra
%           2: ROCK inhibition: factor multiplying param 55, 69 (epsilon and tau - ROCK mediated catalysis)
%           3: MRTF-Circadian coupling inhibition - factor multiplying
%           MRTF-Circadian coupling strength parameters (KeB20,M and KeP20,M)
%           4: YAP overexpression - fold expression of 5SA-YAP (compared to normal YAP expression)
%           5: CytD concentration (in micromolar)
%           6: LATS factor - factor multiplying kNC (rate of YAP phosphorylation)
%           7: blebbistatin - factor multiplying param(46) (rate of stress fiber dependent YAPTAZ
%           dephos) and param(86) (rate of stress fiber-dependent nuclear pore opening)
%           8: cell contact area (in microns squared, control area is 3000)
%           9: lamin A mutation - factor multiplying lamin A phos rate (krl)
%
%     * popVar (optional) - either length of 105 (total # of parameters
%       to vary) or scalar corresponding to population variance

%
% output:
%     * t is the vector of times
%     * y is the vector of state variables
%     * SSVals is the steady-state solution for all state variables in the
%       YAP/TAZ and MRTF mechanotransduction model (see MechanoSS.m)

    if isempty(varargin)
        popVar = 0;
        noiseLevel = [0,0];
    elseif length(varargin)==1
        popVar = varargin{1};
        noiseLevel = [0, 0];
    elseif length(varargin)==2
        popVar = varargin{1};
        noiseLevel = varargin{2};
    end
    
    %Circadian parameters
    tauB = pSol(1);
    pExpB = pSol(2);
    KeB = pSol(3);
    KiB = pSol(4);
    KdBP = pSol(5);
    KdB = pSol(6);
    tauP = pSol(7);
    pExpP = pSol(8);
    KeP = pSol(9);
    KaP = pSol(10);
    KdP = pSol(11);

    magCouple1 = pSol(12); % for YAP/TAZ coupling to BMAL1
    KdCouple1 = pSol(13);
    nCouple1 = pSol(14);
    magCouple2 = pSol(15); % for MRTF coupling to PER/CRY
    KdCouple2 = pSol(16);
    nCouple2 = pSol(17);

    magCouple3 = pSol(20); % for YAP/TAZ coupling to PER/CRY
    KdCouple3 = pSol(21);
    nCouple3 = pSol(22);
    magCouple4 = pSol(23); % for MRTF coupling to BMAL1
    KdCouple4 = pSol(24);
    nCouple4 = pSol(25);

    KdLuc = pSol(29);

    % inhibit MRTF-Circadian coupling if applicable
    magCouple2 = magCouple2 * inhibVec(3);
    magCouple4 = magCouple4 * inhibVec(3);
    
    lags = [tauB, tauP, 1*3600]; % third is nominal delay (only affects time shift of reporter, can be any delay in general)

    SSVals = MechanoSS(stiffnessParam, inhibVec, pSol, popVar);
    YAPTAZnuc_SS = SSVals(15);
    MRTFnuc_SS = SSVals(25);
    mechanoStartEffect = false; % optionally consider effect of YAP/TAZ and MRTF jumping in value after first time lag
    if mechanoStartEffect
        KeB2 = magCouple1 / (1 + (KdCouple1/0.7)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/0.301)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/0.301)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/0.7)^nCouple3);
    else
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3);
    end
    fsolveOptions = optimoptions('fsolve','Display','off');
    unstableSS = fsolve(@(ss) [KeB/(1+(ss(1)/KiB)^pExpB) + KeB2 - KdBP*ss(1)*ss(2) - KdB*ss(1);
                               KeP/(1+(KaP/ss(1))^pExpP) + KeP2 - KdBP*ss(1)*ss(2) - KdP*ss(2)],...
                [0.1,0.1], fsolveOptions);
    unstableSS(3) = 0;
    
    startTime = tic;
    if all(noiseLevel == 0)
        try DDESol = dde23(@ddefun_YAPSS, lags, @history_CircOnly, timeSpan);
            t = DDESol.x';
            y = DDESol.y';
            [~, lucInt] = ode15s(@odeLuc, t, 0, odeset(), t, y(:,1));
            y(:,3) = lucInt';
        catch
            t = (timeSpan(1):3600:timeSpan(end))';
            y = zeros(length(t), 3);
        end
    elseif all(noiseLevel > 0)
        [t, y] = sdde_circadian_solve(noiseLevel(1), noiseLevel(2), noiseLevel(3), lags);
    else
        error('Noise level must be greater than or equal to zero')
    end

    function dydt = ddefun_YAPSS(t,y,Z)
        clockLag1 = Z(1,1);
        clockLag2 = Z(1,2);
        if t < lags(1) && mechanoStartEffect
            YAPTAZnuc_B = 0.7;
            MRTFnuc_B = .301;
        else
            YAPTAZnuc_B = YAPTAZnuc_SS;
            MRTFnuc_B = MRTFnuc_SS;
        end
        if t < lags(2) && mechanoStartEffect
            YAPTAZnuc_P = 0.7;
            MRTFnuc_P = .301;
        else
            YAPTAZnuc_P = YAPTAZnuc_SS;
            MRTFnuc_P = MRTFnuc_SS;
        end
        % clockLag3 = Z(1,3);
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_B)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_B)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_P)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_P)^nCouple3);

        % Rates (BMAL and PER dynamics)
        dydt = [KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdBP*y(1)*y(2) - KdB*y(1) + KeB2;
                KeP*(1./(1+(KaP/clockLag2).^pExpP)) - KdBP*y(1)*y(2) - KdP*y(2) + KeP2];
                % KeP*(1./(1+(KaP/clockLag3).^pExpP)) - KdLuc*y(3) + KeP2];
        curEval = toc(startTime);
        if curEval > 10 % return error if this instance has taken longer than 10 s to run
            error('Solver took too long')
        end
    end

    function s = history_CircOnly(~)
        s = [5*unstableSS(1); 0.2*unstableSS(2)];%; unstableSS(3)];
    end

    function dluc = odeLuc(t, luc, tStored, BStored)
        BCur = interp1(tStored, BStored, t);
        dluc = KeP*(1./(1+(KaP/BCur).^pExpP)) - KdLuc*luc + KeP2;
    end

    function [t,y] = sdde_circadian_solve(eps_B, eps_P, eps_L, lags)
        % solve stochastic version of this system of ddes
        % not used in the current paper

        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3);
        dt = 2; % ~2 s time step
        delayFactor1 = round(lags(1)/dt);
        delayFactor2 = round(lags(2)/dt);
        t = (0:dt:timeSpan(2))';
        y = zeros(length(t),3);
        y(1,1:2) = history_CircOnly(0);
        y(1,3) = 0;
        dy = zeros(1,3);
        eta_B = 0;
        eta_P = 0;
        eta_L = 0;
        rand_stored = sqrt(dt)*randn(length(t)-1, 3);
        for i = 2:length(t)
            eta_B = eta_B - dt*eta_B + eps_B*rand_stored(i-1,1) + 0.5*eps_B^2*(rand_stored(i-1,1)^2-dt);
            eta_P = eta_P - dt*eta_P + eps_P*rand_stored(i-1,2) + 0.5*eps_P^2*(rand_stored(i-1,2)^2-dt);
            eta_L = eta_L - dt*eta_L + eps_L*rand_stored(i-1,3) + 0.5*eps_L^2*(rand_stored(i-1,3)^2-dt);
            % Rates (BMAL and PER dynamics)
            B = y(i-1,1);
            P = y(i-1,2);
            L = y(i-1,3);
            if t(i-1) < lags(1)
                BLag1 = y(1,1);
            else
                BLag1 = y(i-1-delayFactor1,1);
            end
            if t(i-1) < lags(2)
                BLag2 = y(1,1);
            else
                BLag2 = y(i-1-delayFactor2,1);
            end
            dy(1) = KeB*(1./(1+(BLag1/KiB).^pExpB)) - KdBP*B*P - KdB*B + KeB2 + eta_B;
            dy(2) = KeP*(BLag2^pExpP./(BLag2^pExpP+KaP^pExpP)) - KdBP*B*P - KdP*P + KeP2 + eta_P;
            dy(3) = KeP*(1./(1+(KaP/B).^pExpP)) - KdLuc*L + KeP2 + eta_L;
            y(i,:) = y(i-1,:) + dy*dt;
            if any(y(i,:) < 0)
                yCur = y(i,:);
                yCur(yCur < 0) = 0;
                y(i,:) = yCur;
            end
        end
    end
end