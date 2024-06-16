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
% INPUTS:
%     * pSol: parameter solution vector - length 45, contains all calibrated
%           values from mechano-Circadian model (see Tables S1 and S4 in paper).
%           In order, the parameters are:
%            {'tauB','nB','KeB0','KiB','KdBP','KdB',...
%             'tauP','nP0','KeP0','KaP','KdP',...
%             'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
%             'ROCKInhibSens', 'StiffThresh',...
%             'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
%             'LatBSens','JasMag','JasSens','KdLuc',...
%             'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim','Kcap',...
%             'KeP1', 'KiP', 'KdR', 'tauR',...
%             'KeR2,Y', 'KYR', 'KeR2,M', 'KMR', 'nP1', 'nYR', 'nMR'};
%       * inhibVec: vector of inhibition parameters (length=9 or 10)
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
%           10: (optional) lamin A mutation - factor multiplying NPC opening rate (KfNPC)
%
%       Additionally, we have the following optional arguments (set to
%       defaults if not provided)
%       * popVar: either length of 105 (total # of parameters
%           to vary) or scalar corresponding to population variance
%           (default is 0)
%       * noiseLevel: noise parameters associated with SDDE implementation
%           [noiseB, noiseP, noiseR], default is [0,0,0]
%       * Ke3: baseline expression independent of YAP/TAZ and MRTF
%           [Ke3_B, Ke3_P, Ke3_R], default is [0,0,0]

%
% output:
%     * t is the vector of times
%     * y is the vector of state variables
%     * SSVals is the steady-state solution for all state variables in the
%       YAP/TAZ and MRTF mechanotransduction model (see MechanoSS.m)

    if isempty(varargin)
        popVar = 0;
        noiseLevel = [0,0];
        Ke3 = [0,0,0];
    elseif length(varargin)==1
        popVar = varargin{1};
        noiseLevel = [0, 0];
        Ke3 = [0,0,0];
    elseif length(varargin)==2
        popVar = varargin{1};
        noiseLevel = varargin{2};
        Ke3 = [0,0,0];
    elseif length(varargin)==3
        popVar = varargin{1};
        noiseLevel = varargin{2};
        Ke3 = varargin{3};
    end
    
    %Circadian parameters
    tauB = pSol(1);
    pExpB = pSol(2);
    KeB = pSol(3);
    KiB = pSol(4);
    % KdBP = pSol(5);
    KdB = pSol(6);
    tauP = pSol(7);
    pExpP0 = pSol(8);
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

    
    KeP0_self = pSol(35);
    KiP = pSol(36);
    KdR = pSol(37);
    tauR = pSol(38);
    pExpP1 = pSol(43);
    magCouple5 = pSol(39); % for YAP/TAZ coupling to REV-ERBalpha
    KdCouple5 = pSol(40);
    nCouple5 = pSol(44);
    magCouple6 = pSol(41); % for MRTF coupling to REV-ERBalpha
    KdCouple6 = pSol(42);
    nCouple6 = pSol(45);

    % inhibit MRTF-Circadian coupling if applicable
    magCouple2 = magCouple2 * inhibVec(3);
    magCouple4 = magCouple4 * inhibVec(3);
    
    lags = [tauB, tauP, tauR, 1*3600]; % 4th lag is arbitrary, just determines phase of reporter

    SSVals = MechanoSS(stiffnessParam, inhibVec, pSol, popVar);
    YAPTAZnuc_SS = SSVals(15);
    MRTFnuc_SS = SSVals(25);
    mechanoStartEffect = false; % optionally consider effect of YAP/TAZ and MRTF jumping in value after first time lag
    if mechanoStartEffect
        KeB2 = magCouple1 / (1 + (KdCouple1/0.7)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/0.301)^nCouple4) + Ke3(1);
        KeP2 = magCouple2 / (1 + (KdCouple2/0.301)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/0.7)^nCouple3) + Ke3(2);
        KeR2 = magCouple5 / (1 + (KdCouple5/0.7)^nCouple5) +...
               magCouple6 / (1 + (KdCouple6/0.301)^nCouple6) + Ke3(3);
    else
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4) + Ke3(1);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3) + Ke3(2);
        KeR2 = magCouple5 / (1 + (KdCouple5/YAPTAZnuc_SS)^nCouple5) +...
               magCouple6 / (1 + (KdCouple6/MRTFnuc_SS)^nCouple6) + Ke3(3);
    end

    fsolveOptions = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12);
    initGuess = [1, 1, 1];
    unstableSS = fsolve(@(ss) [KeB/(1+(ss(3)/KiB)^pExpB) + KeB2 - KdB*ss(1);
                               KeP/(1+(KaP/ss(1))^pExpP0) + KeP0_self*(1./(1+(ss(2)/KiP).^pExpP1)) + KeP2 - KdP*ss(2);
                               KeP/(1+(KaP/ss(1))^pExpP0) + KeP0_self*(1./(1+(ss(2)/KiP).^pExpP1)) + KeR2 - KdR*ss(3)],...
                initGuess, fsolveOptions);

    if any(unstableSS < 0)
        fprintf("Warning: negative initial condition being set to zero\n")
        unstableSS(unstableSS < 0) = 0;
    end
    startTime = tic;
    if all(noiseLevel == 0)
        try DDESol = dde23(@ddefun_YAPSS, lags, @history_CircOnly, timeSpan);
            t = DDESol.x';
            y = DDESol.y';
            [~, lucInt] = ode15s(@odeLuc, t, 0, odeset(), t, y);
            y(:,4) = lucInt';
        catch
            t = (timeSpan(1):3600:timeSpan(end))';
            y = zeros(length(t), 4);
        end
    else
        error('Stochastic version not currently supported')
    end

    function dydt = ddefun_YAPSS(t,y,Z)
        clockLag1 = Z(3,1);%(1,1);
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
        if t < lags(3) && mechanoStartEffect
            YAPTAZnuc_R = 0.7;
            MRTFnuc_R = .301;
        else
            YAPTAZnuc_R = YAPTAZnuc_SS;
            MRTFnuc_R = MRTFnuc_SS;
        end
        clockLag3 = Z(2,2);
        clockLag4 = Z(1,3);
        clockLag5 = Z(2,3);
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_B)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_B)^nCouple4) + Ke3(1);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_P)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_P)^nCouple3) + Ke3(2);
        KeR2 = magCouple5 / (1 + (KdCouple5/YAPTAZnuc_R)^nCouple5) +...
               magCouple6 / (1 + (KdCouple6/MRTFnuc_R)^nCouple6) + Ke3(3);

        % Rates (BMAL and PER dynamics)
        dydt = [KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdB*y(1) + KeB2;
                KeP*(1./(1+(KaP/clockLag2).^pExpP0)) + KeP0_self*(1./(1+(clockLag3/KiP).^pExpP1)) ...
                - KdP*y(2) + KeP2;
                KeP*(1./(1+(KaP/clockLag4).^pExpP0)) + KeP0_self*(1./(1+(clockLag5/KiP).^pExpP1)) ...
                - KdR*y(3) + KeR2];
        curEval = toc(startTime);
        if curEval > 10 % return error if this instance has taken longer than 10 s to run
            error('Solver took too long')
        end
    end

    function s = history_CircOnly(~)
        s = [2*unstableSS(1); 0.2*unstableSS(2); unstableSS(3)];
    end

    function dluc = odeLuc(t, luc, tStored, yVals)
        if t > lags(4)
            BCur = interp1(tStored, yVals(:,1), t-lags(4));
            PCur = interp1(tStored, yVals(:,2), t-lags(4));
        else
            BCur = yVals(1,1);
            PCur = yVals(1,2);
        end
        KeP_self = KeP0_self*(1./(1+(PCur/KiP).^pExpP1));
        dluc = KeP*(1./(1+(KaP/BCur).^pExpP0)) + KeP_self - KdLuc*luc + KeP2;
    end

end