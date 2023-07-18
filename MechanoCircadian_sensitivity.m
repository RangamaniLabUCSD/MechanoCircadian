clearvars
rng(100,'twister')
uqlab
ModelOpts.mFile = 'uq_MechanoCircadian';
myModel = uq_createModel(ModelOpts);
paramNames = {'tau','nB','KeB0','KiB','KdBP','KdB',...
                   'nP','KeP0','KaP','KdP',...
              'KeB2,Y','KYB','KeP2,M','KMP',...
              'KeP2,Y','KYP','KeB2,M','KMB'};
% tauB1 = p(1);     1
% pExpB = p(2);     2
% KeB = p(3);       3
% KiB = p(4);       4
% KdBP = p(5);      5
% KdB = p(6);       6
% tauB2 = p(7);     -
% pExpP = p(8);     7
% KeP = p(9);       8
% KaP = p(10);      9
% KdP = p(11);      10
% magCouple1 = p(12);11
% KdCouple1 = p(13);12
% nCouple1 = p(14);-
% magCouple2 = p(15);13
% KdCouple2 = p(16);14
% nCouple2 = p(17);-
% % p(18) is CytD sens.
% C = p(19);%stiffness shift
% 20: magCouple3    15
% 21: KdCouple3     16
% 22: nCouple3      
% 23: magCouple4    17
% 24: KdCouple4     18
% 25: nCouple4
% 26: LatB const
% 27: Jasp mag const
% 28: Jasp sens const
numParam = 18;
for i = 1:numParam
    InputOpts.Marginals(i).Type = 'Lognormal';
    InputOpts.Marginals(i).Name = paramNames{i};
    InputOpts.Marginals(i).Parameters = [0 log(1.5)];
end
myInput = uq_createInput(InputOpts);
SobolOpts.Type = 'uq_sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 4096;
SobolOpts.Bootstrap.Replications = 100;
SobolOpts.Bootstrap.Alpha = 0.05;
SobolAnalysis = uq_createAnalysis(SobolOpts);

%% display
uq_visualizeSobolIndices(SobolAnalysis.Results, 1:3)