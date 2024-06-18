function [periodTest, amplTest, tOut, yOut, rawOutput, oscDecayRate] = ...
    conditionToOutputs(p, stiffness, inhibVec, varargin)
    % function taking in conditions for a MechanoCircadian test and
    % outputting the oscillation period, amplitude, and decay rate, along
    % with the oscillation dynamics.

    % INPUTS:
    %   * p: parameter solution vector - length 45, contains all calibrated
%           values from mechano-Circadian model (see Tables S1 and S4 in paper).
%           In order, the parameters are:
%           {'tauB','nB','KeB0','KiB','KdBP','KdB',...
%             'tauP','nP0','KeP0','KaP','KdP',...
%             'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
%             'ROCKInhibSens', 'StiffThresh',...
%             'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
%             'LatBSens','JasMag','JasSens','KdLuc',...
%             'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim','Kcap',...
%             'KeP1', 'KiP', 'KdR', 'tauR',...
%             'KeR2,Y', 'KYR', 'KeR2,M', 'KMR', 'nP1', 'nYR', 'nMR'};
%       * stiffness: substrate stiffness in kPa
%       * inhibVec: vector of inhibition parameters (length=9 or 10)
%           set to default if empty vector []
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
%       Additionally, we have the following optional arguments (set to
%       defaults if not provided)
%       * maxTime: max time for the simulation in s (default is 2400*24*7 -
%           1 week)
%       * popVar: either length of 105 (total # of parameters
%           to vary) or scalar corresponding to population variance
%           (default is 0)
%       * noiseLevel: noise parameters associated with SDDE implementation
%           [noiseB, noiseP, noiseR], default is [0,0,0]
%       * Ke3: baseline expression independent of YAP/TAZ and MRTF
%           [Ke3_B, Ke3_P, Ke3_R], default is [0,0,0]

%   OUTPUTS:
%       * periodTest: oscillation period in s for each of the 4 MechanoCircadian variables, 
%         BMAL1 (B), PER/CRY (P), REV-ERBalpha (R), and the PER2 luciferase reporter (L)
%         (determined using circOscAnalysis)
%       * amplTest: oscillation amplitude for each of the 4 MechanoCircadian variables 
%         (determined using circOscAnalysis)
%       * tOut: relative time output (aligned such that second peak in
%       luciferase reporter occurs at time t=0), time interval is 3600 s
%       * yOut: output for each state variable at times tOut (4 columns) 
%       * rawOutput: cell containing raw outputs from MechanoCircadianModel
%       - {t, y, ySS}
%       * oscDecayRate: oscillation period in s for each of the 4 MechanoCircadian variables 
%         (determined using circOscAnalysis)

    % handle different numbers of optional args
    if isempty(varargin)
        maxTime = 3600*24*7;
        popVar = 0;
        noiseLevel = [0,0];
        Ke3 = [0,0,0];
    elseif length(varargin)==1
        maxTime = varargin{1};
        popVar = 0;
        noiseLevel = [0,0];
        Ke3 = [0,0,0];
    elseif length(varargin)==2
        maxTime = varargin{1};
        popVar = varargin{2};
        noiseLevel = [0,0];
        Ke3 = [0,0,0];
    elseif length(varargin)==3
        maxTime = varargin{1};
        popVar = varargin{2};
        noiseLevel = varargin{3};
        Ke3 = [0,0,0];
    elseif length(varargin)==4
        maxTime = varargin{1};
        popVar = varargin{2};
        noiseLevel = varargin{3};
        Ke3 = varargin{4};
    end
    % run MechanoCircadianMOdel for current conditions
    [t,y,ySS] = MechanoCircadianModel([0 maxTime], stiffness, p, inhibVec, popVar, noiseLevel, Ke3);
    % Return default vector if any nonreal or negative solutions
    if any(~isreal(y(:))) || any(y(:)<0)
        periodTest = zeros(1,size(y,2));
        amplTest = 1e-6*ones(1,size(y,2));
        tOut = (0:3600:maxTime/2)';
        yOut = zeros([length(tOut), 4]);
        rawOutput = {t, y, ySS};
        oscDecayRate = -1*ones(1,size(y,2));
        return
    end
    periodTest = zeros(1, size(y,2));
    amplTest = zeros(1, size(y,2));
    oscDecayRate = zeros(1, size(y,2));
    locs = cell(1, size(y,2));
    for i = 1:size(y,2)
        oscDynamics = y(:,i);
        [periodTest(i), amplTest(i), oscDecayRate(i), locs{i}] = circOscAnalysis(t, oscDynamics);
    end
    % interpolate dynamics with t=0 defined at the second peak in the
    % luciferase reporter (L), or at the maximum value of L if there are 2
    % or less peaks overall
    if length(locs{4})>2
        maxTimeOut = min([maxTime*0.9, t(end)-t(locs{4}(2))]);
        tOut = (0:3600:maxTimeOut)';
        yOut = interp1(t-t(locs{4}(2)), y, tOut);
    else
        idx = find(t < 0.1*maxTime, 1, 'last');
        [~,maxIdx] = max(y(1:idx,4));
        tOut = (0:3600:maxTime*0.9)';
        yOut = interp1(t-t(maxIdx), y, tOut);
    end
    rawOutput = {t, y, ySS};
end