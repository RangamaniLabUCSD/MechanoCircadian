function [period, ampl, decayRate, pkLocs] = circOscAnalysis(t, oscDynamics)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    maxTime = max(t);
    [pks,locs] = findpeaks(oscDynamics,'MinPeakProminence',.001);
    if length(pks) > 2
        periodEst = mean(diff(t(locs(2:end)))); % rough estimate of period
        if periodEst < 2*3600 % should not be less than 2 hours!
            periodEst = 24*3600;
        end
    else
        period = t(end);
        ampl = max(oscDynamics)-min(oscDynamics);
        decayRate = -1;
        pkLocs = [];
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
    pkLocs = zeros([length(passThroughMeanIdx)-1,1]);
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
            pkLocs(i) = idx + curLoc - 1;
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
            smallOscIdx = find(amplitudeStored > 0.1*amplitudeStored(1),1,'last');
            smallOscIdx = max([2, smallOscIdx]);
            period = mean(diff(tZeros(2:smallOscIdx+1)));
            ampl = mean(amplitudeStored(1:smallOscIdx));
            if amplitudeStored(1) > 0
                decayRate = 3600*mean(diff(amplitudeStored(1:smallOscIdx))./...
                    diff(tInterp(pkLocs(2:smallOscIdx+1)))')/amplitudeStored(1);
            else
                decayRate = -1;
            end
        else
            if amplitudeStored(1) > 0
                decayRate = 3600*mean(diff(amplitudeStored)./...
                    diff(tInterp(pkLocs(2:end)))')/amplitudeStored(1);
            else
                decayRate = -1;
            end
            period = mean(diff(tZeros(2:end)));
            ampl = mean(amplitudeStored);
        end
    else
        decayRate = -1;
        period = t(end);
        ampl = max(oscDynamics)-min(oscDynamics);
    end
    % convert to pkLocs in terms of original time vector
    for i = 1:length(pkLocs)
        [~,pkLocs(i)] = min(abs(tInterp(pkLocs(i))-t));
    end
end