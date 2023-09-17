function [period, ampl, decayRate, pkLocs] = circOscAnalysis(t, oscDynamics)
    % function for computing oscillation period, amplitude, and decay rate
    % for MechanoCircadian model

    % INPUTS:
%       * t: time vector
%       * oscDynamics: vector storing values for variable of interest

    % OUTPUTS:
%       * period: oscillation period
%       * ampl: oscillation amplitude
%       * decayRate: oscillation decay rate (essentially, the rate of
%       change in amplitude)
%       * pkLocs: locations of peaks, given in terms of indices

    maxTime = max(t);
    [pks,locs] = findpeaks(oscDynamics,'MinPeakProminence',mean(oscDynamics)*1e-4);
    if length(pks) > 2
        periodEst = mean(diff(t(locs(2:end)))); % rough estimate of period
        if periodEst < 2*3600 % should not be less than 2 hours!
            periodEst = 24*3600;
        end
    else
        period = 0;%t(end);
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
    yInterp = smooth(yInterp,round(0.25*periodEst/dt));
    % [~,locsSmooth] = findpeaks(smooth(yInterp,101),'MinPeakProminence',.001);
    meanDyn = 0;%mean(yInterp(tInterp>=maxTime/2));
    passThroughMeanIdx = find((yInterp(1:end-1)<=meanDyn & yInterp(2:end)>meanDyn));
    passThroughMeanIdx = passThroughMeanIdx(passThroughMeanIdx > filtWindow/2 ...
        & passThroughMeanIdx < (length(yInterp)-filtWindow/2));

    tZeros = zeros(size(passThroughMeanIdx));
    pks = zeros([length(passThroughMeanIdx)-1,1]);
    pkLocs = zeros([length(passThroughMeanIdx)-1,1]);
    troughs = zeros([length(passThroughMeanIdx)-1,1]);
    troughLocs = zeros([length(passThroughMeanIdx)-1,1]);
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
            [~,curTroughLoc] = min(yInterp(idx:nextIdx));
            troughLocs(i) = idx + curTroughLoc - 1;
        end
    end

    if length(pks)>2 && length(troughs)>2
        lastExtremum = min([length(pks),length(troughs)]); 
        amplitudeStored = pks(1:lastExtremum) - (troughs(1:lastExtremum));
        if amplitudeStored(end) < 0.5*amplitudeStored(1) % detect special case of decaying oscillations
            smallOscIdx = find(amplitudeStored > 0.5*amplitudeStored(1),1,'last');
            smallOscIdx = max([2, smallOscIdx]);
            period = mean(diff(tZeros(1:smallOscIdx+1)));
            % gapEst = [diff(tInterp(pkLocs(2:smallOscIdx+1))), diff(tInterp(troughLocs(2:smallOscIdx+1)))];
            % period = mean(gapEst);
            ampl = mean(amplitudeStored(1:smallOscIdx));
            if amplitudeStored(1) > 0
                decayRate = 3600*mean(diff(amplitudeStored(1:smallOscIdx))./...
                    diff(tInterp(pkLocs(1:smallOscIdx)))')/amplitudeStored(1);
            else
                decayRate = -1;
            end
        else
            if amplitudeStored(1) > 0
                decayRate = 3600*mean(diff(amplitudeStored)./...
                    diff(tInterp(pkLocs(1:end)))')/amplitudeStored(1);
            else
                decayRate = -1;
            end
            period = mean(diff(tZeros(1:end)));
            ampl = mean(amplitudeStored);
        end
    else
        decayRate = -1;
        period = 0;%t(end);
        ampl = max(oscDynamics)-min(oscDynamics);
    end
    % convert to pkLocs in terms of original time vector
    for i = 1:length(pkLocs)
        [~,pkLocs(i)] = min(abs(tInterp(pkLocs(i))-t));
    end
end