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
    qoi = zeros(numSamples,3);
    maxTime = 480*3600;
    % p0 = [8*3600; 2.5; 1/3600; .04; 0.5/3600; 0.4/3600; 8*3600; 2.5; 1/3600; .5; .4/3600;...
    %       0.05/3600; 1; 2; .05/3600; 1; 2; 1; 3.25];
    p0 = [29003.5906784085	2.50000000000000	0.00245158891914641	0.364964919887603	...
        0.000159446594124508	4.35084376960400e-05	28800	2.50000000000000	0.000227497666625881	...
        1.30746956455001	0.000559475102598622	0.000138033797382737	3.63884412637135	2	...
        0.000134584220775489	7.62727798452122	2	99.8862115161598	14.4863689973732	...
        1.78509933795973e-05	5.85622066745939	2	8.81762812384363e-05	8.33653139841384	...
        2	1.94985574208113	13.7231734117393	0.999997028946840]';
    activeIdx = [1:6, 8:13, 15,16,20,21,23,24];
    parfor i = 1:numSamples
        pCur = p0;
        pCur(activeIdx) = p(i,:).*p0(activeIdx)';
        inhibVec = [1, 1, 1, 1, 0];
        [curPeriod, curAmpl,~,~,~, oscDecayRate] = conditionToOutputs(pCur,1e7,inhibVec,maxTime);
        qoi(i,:) = [curPeriod, curAmpl, oscDecayRate];
    end
end