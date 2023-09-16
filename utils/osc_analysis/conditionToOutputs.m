function [periodTest, amplTest, tOut, yOut, rawOutput, oscDecayRate] = conditionToOutputs(p, stiffness, inhibVec, varargin)
    if isempty(varargin)
        maxTime = 3600*24*7;
        popVar = 0;
        noiseLevel = [0,0];
    elseif length(varargin)==1
        maxTime = varargin{1};
        popVar = 0;
        noiseLevel = [0,0];
    elseif length(varargin)==2
        maxTime = varargin{1};
        popVar = varargin{2};
        noiseLevel = [0,0];
    elseif length(varargin)==3
        maxTime = varargin{1};
        popVar = varargin{2};
        noiseLevel = varargin{3};
    end
    [t,y,ySS] = MechanoCircadianModel([0 maxTime], [stiffness,inf,0], p, inhibVec, popVar, noiseLevel);
    if any(~isreal(y(:))) || any(y(:)<0)
        periodTest = zeros(1,size(y,2));
        amplTest = 1e-6*ones(1,size(y,2));
        tOut = (0:3600:maxTime/2)';
        yOut = zeros([length(tOut), 3]);
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
    if length(locs{3})>2
        maxTimeOut = min([maxTime*0.9, t(end)-t(locs{3}(2))]);
        tOut = (0:3600:maxTimeOut)';
        yOut = interp1(t-t(locs{3}(2)), y, tOut);
    else
        tOut = (0:3600:maxTime*0.9)';
        yOut = interp1(t, y, tOut);
    end
    rawOutput = {t, y, ySS};
end