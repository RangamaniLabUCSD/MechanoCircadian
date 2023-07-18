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
    if any(~isreal(y))
        periodTest = maxTime;
        amplTest = 1e-6;
        tOut = (0:3600:maxTime/2)';
        yOut = zeros([length(tOut), 3]);
        rawOutput = {t, y, ySS};
        oscDecayRate = -1;
        return
    end
    % YAPTAZnuc_SS = ySS(15);
    % MRTFnuc_SS = ySS(25);
    % KeB2 = p(12) / (1 + (p(13)/YAPTAZnuc_SS)^2) +...
    %        p(23) / (1 + (p(24)/MRTFnuc_SS)^2);
    % KeP2 = p(15) / (1 + (p(16)/MRTFnuc_SS)^2) +...
    %        p(20) / (1 + (p(21)/YAPTAZnuc_SS)^2);
    oscDynamics = y(:,3);%3600*p(9)./(1 + (p(10)./y(:,1)).^p(8)) + KeB2 + KeP2;
    % y(:,3) = oscDynamics;
    [pks,locs] = findpeaks(oscDynamics,'MinPeakProminence',.001);
    if length(pks) > 2
        periodEst = mean(diff(t(locs(2:end)))); % rough estimate of period
        if periodEst < 2*3600 % should not be less than 2 hours!
            periodEst = 24*3600;
        end
    else
        periodTest = t(end);
        amplTest = max(oscDynamics)-min(oscDynamics);
        oscDecayRate = -1;
        tOut = (0:3600:maxTime*0.9)';
        yOut = interp1(t, y, tOut);
        rawOutput = {t, y, ySS};
        return
    end
    dt = periodEst/500;
    tInterp = 0:dt:maxTime; % sampling about every 1/500 of the period
    yInterp = interp1(t,smooth(oscDynamics),tInterp);
    filtWindow = round(periodEst/dt);
    yInterp = yInterp' - smooth(yInterp,filtWindow);
    % [~,locsSmooth] = findpeaks(smooth(yInterp,101),'MinPeakProminence',.001);
    meanDyn = 0;%mean(yInterp(tInterp>=maxTime/2));
    passThroughMeanIdx = find((yInterp(1:end-1)<=meanDyn & yInterp(2:end)>meanDyn));
    passThroughMeanIdx = passThroughMeanIdx(passThroughMeanIdx > filtWindow/2 ...
        & passThroughMeanIdx < (length(yInterp)-filtWindow/2));

    tZeros = zeros(size(passThroughMeanIdx));
    pks = zeros([length(passThroughMeanIdx)-1,1]);
    locs = zeros([length(passThroughMeanIdx)-1,1]);
    troughs = zeros(size(passThroughMeanIdx));
    for i = 1:length(passThroughMeanIdx)
        idx = passThroughMeanIdx(i);
        tZeros(i) = tInterp(idx) + (tInterp(idx+1)-tInterp(idx)) * (meanDyn-yInterp(idx))/(yInterp(idx+1)-yInterp(idx));
        if i < length(passThroughMeanIdx)
            nextIdx = passThroughMeanIdx(i+1);
            oscDynRange = t >= tInterp(idx) & t < tInterp(nextIdx);
            if sum(oscDynRange)>0
                pks(i) = max(oscDynamics(oscDynRange));
                troughs(i) = min(oscDynamics(oscDynRange));
            else
                curVals = interp1(t,oscDynamics,[tInterp(idx), tInterp(nextIdx)]);
                pks(i) = max(curVals);
                troughs(i) = min(curVals);
            end
            [~,curLoc] = max(yInterp(idx:nextIdx));
            locs(i) = idx + curLoc - 1;
        end
    end

    % fft to find peak period?
%     tInterp = 0:3600:maxTime;
%     tInterp = tInterp(tInterp > t(locs(2)));
%     yInterp = interp1(t,oscDynamics,tInterp);
%     yInterp = yInterp - mean(yInterp); % remove zero-order
%     Fs = 1/3600;
%     N = length(tInterp);
%     Y_fft = fft(yInterp);
%     Y_fft = Y_fft(1:floor(N/2+1));
%     curPower = (1/(Fs*N)) * abs(Y_fft).^2;
%     curPower(2:end-1) = 2*curPower(2:end-1);
%     freq = 0:Fs/N:Fs/2;
%     [~,maxFreqIdx] = max(curPower);
%     periodTestFFT = 1/freq(maxFreqIdx);
    if length(pks)>2 && length(troughs)>2
        lastExtremum = min([length(pks),length(troughs)]); 
        amplitudeStored = pks(2:lastExtremum) - (troughs(2:lastExtremum));
        if amplitudeStored(end) < 0.5*amplitudeStored(1) % if it is decaying, mark as nonosc.
            oscDecayRate = 3600*mean(diff(amplitudeStored)./...
                diff(tInterp(locs(2:end)))')/amplitudeStored(1);
            smallOscIdx = find(amplitudeStored > 0.1*amplitudeStored(1),1,'last');
            periodTest = mean(diff(tZeros(2:smallOscIdx+1)));
            amplTest = mean(amplitudeStored(1:smallOscIdx));
        else
            oscDecayRate = 3600*mean(diff(amplitudeStored)./...
                diff(tInterp(locs(2:end)))')/amplitudeStored(1);
            periodTest = mean(diff(tZeros(2:end)));
            amplTest = mean(amplitudeStored);
        end
        maxTimeOut = min([maxTime*0.9, t(end)-tInterp(locs(2))]);
        tOut = (0:3600:maxTimeOut)';
        yOut = interp1(t-tInterp(locs(2)), y, tOut);
    else
        oscDecayRate = -1;
        periodTest = t(end);
        amplTest = max(oscDynamics)-min(oscDynamics);
        tOut = (0:3600:maxTime*0.9)';
        yOut = interp1(t, y, tOut);
    end
    rawOutput = {t, y, ySS};
end