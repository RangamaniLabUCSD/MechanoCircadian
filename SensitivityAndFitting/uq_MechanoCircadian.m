function qoi = uq_MechanoCircadian(p)
    % wrapper function for sensitivity analysis using UQLab
    % INPUTS:
    %   p: matrix of parameters, where each row represents a sample and
    %       each column represents a parameter (e.g. for the first sample,
    %       third parameter, use p(1,3)
    %
    % OUTPUTS:
    %   qoi: matrix of quantities of interest (QoIs) for each sample. In
    %   this case, there are 36 quantities of interest - 4 stiffnesses are
    %   tested (1e7, 100, 10, 1 kPa), and for each stiffness, we record the
    %   period, amplitude, and decay rate for oscillations associated with
    %   each Circadian variable (BMAL1, PER/CRY, and the luciferase
    %   reporter)
    %
    % Note that this function currently uses a parallel for loop with 12
    % workers to speed up sensitivity analysis. 

    numSamples = size(p,1);
    qoi = zeros(numSamples,48);
    maxTime = 480*3600;
    p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
        0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
        100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
    activeIdx = [1:4, 6:13, 15,16,19,20,21,23,24,29,30,31,32,35:43];
    stiffnessVals = [1e7, 100, 10, 1];
    parfor i = 1:numSamples
        pCur = p0;
        pCur(activeIdx) = p(i,:).*p0(activeIdx)';
        inhibVec = [1, 1, 1, 0, 0, 1];
        inhibVec(7:9) = [1,3000,1];
        qoiCur = zeros(1,12*length(stiffnessVals));
        for j = 1:length(stiffnessVals)
            [curPeriod, curAmpl,~,~,~, oscDecayRate] = conditionToOutputs(pCur,stiffnessVals(j),inhibVec,maxTime);
            qoiCur((j-1)*12+1:j*12) = [curPeriod, curAmpl, oscDecayRate];
        end
        qoi(i,:) = qoiCur;
        if ~mod(i,100)
            fprintf('Finished step %d of %d\n', i, numSamples)
        end
    end
end