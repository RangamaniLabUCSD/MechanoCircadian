clearvars
rng(100,'twister')
uqlab
ModelOpts.mFile = 'uq_MechanoCircadian';
myModel = uq_createModel(ModelOpts);
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
% % p(18) is old, unused
% C = p(19);%stiffness shift        15
numParam = 15;
for i = 1:numParam
    InputOpts.Marginals(i).Type = 'Lognormal';
    InputOpts.Marginals(i).Parameters = [0 log(1.25)];
end
myInput = uq_createInput(InputOpts);
SobolOpts.Type = 'uq_sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 1024;
SobolOpts.Bootstrap.Replications = 100;
SobolOpts.Bootstrap.Alpha = 0.05;
SobolAnalysis = uq_createAnalysis(SobolOpts);

