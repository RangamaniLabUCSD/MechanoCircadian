function qoi = uq_MechanoCircadian(p)
    % wrapper function for sensitivity analysis using UQLab
    % INPUTS:
    %   p: matrix of parameters, where each row represents a sample and
    %       each column represents a parameter (e.g. for the first sample,
    %       third parameter, use p(1,3)
    %
    % OUTPUTS:
    %   qoi: matrix of quantities of interest (QoIs) for each sample

    numSamples = size(p,1);
    qoi = zeros(numSamples,2);
    maxTime = 500*3600;
    % p0 = [8*3600; 2.5; 1/3600; .04; 0.5/3600; 0.4/3600; 8*3600; 2.5; 1/3600; .5; .4/3600;...
    %       0.05/3600; 1; 2; .05/3600; 1; 2; 1; 3.25];
    p0 = [30112.1128026662	2.49361433104013	0.00110195319176539	0.0362906811897288	...
        0.000138718472336416	4.01740850738952e-05	28800	1.13993216317370	0.000138008691550797	...
        0.316725987150128	0.000225606126486809	4.35158691215436e-05	3.94537576531232	1.86821956003706    ...
    	1.19993905179384e-05	3.95587920932383	1.60484588806550	2.37073084328230	10.2086387270683]';
    activeIdx = [1:6, 8:13, 15,16,19];
    for i = 1:numSamples
        pCur = p0;
        pCur(activeIdx) = p(i,:).*p0(activeIdx)';
        inhibVec = [1, 1, 1, 1, 0];
        [t,y] = MechanoCircadianModel([0 maxTime], [1e5,inf,0], pCur, inhibVec, 0);
        if any(~isreal(y))
            qoi(i,1) = maxTime;
            qoi(i,2) = 0;
        else
            [pks,locs] = findpeaks(y(:,2),'MinPeakProminence',.01);
            [troughs,~] = findpeaks(-y(:,2),'MinPeakProminence',.01);
            if length(pks)<=2 || length(troughs)<=2
                qoi(i,1) = t(end);
                qoi(i,2) = 0;
            else
                lastExtremum = min([length(pks),length(troughs)]); 
                amplitudeStored = pks(2:lastExtremum) - (-troughs(2:lastExtremum));
                qoi(i,2) = mean(amplitudeStored);
                if amplitudeStored(end) < 0.8*amplitudeStored(1) % cannot be decaying rapidly (not sustained osc)
                    qoi(i,1) = t(end);
                else
                    qoi(i,1) = mean(diff(t(locs(2:end))));
                end
            end
        end
    end
end