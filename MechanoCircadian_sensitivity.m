clearvars
rng(100,'twister')
uqlab
ModelOpts.mFile = 'uq_MechanoCircadian';
myModel = uq_createModel(ModelOpts);
paramNames = {'tauB','nB','KeB0','KiB','KdBP','KdB',...
              'tauP','nP','KeP0','KaP','KdP',...
              'KeB2,Y','KYB','nYB','KeP2,M','KMP','nMP'...
              'Kcap', 'StiffThresh',...
              'KeP2,Y','KYP','nYP','KeB2,M','KMB','nMB',...
              'LatBSens','JasMag','JasSens','KdLuc',...
              'MRTFRelease','KinsoloMRTF','Kin2MRTF','Kdim'};
% tauB = p(1);     1
% pExpB = p(2);     
% KeB = p(3);       2
% KiB = p(4);       3
% KdBP = p(5);      4
% KdB = p(6);       5
% tauP = p(7);     6
% pExpP = p(8);     
% KeP = p(9);       7
% KaP = p(10);      8
% KdP = p(11);      9
% magCouple1 = p(12);10
% KdCouple1 = p(13);11
% nCouple1 = p(14);-
% magCouple2 = p(15);12
% KdCouple2 = p(16);13
% nCouple2 = p(17);-
% % p(18) is CytD sens.
% C = p(19);%stiffness shift 14
% 20: magCouple3    15
% 21: KdCouple3     16
% 22: nCouple3      
% 23: magCouple4    17
% 24: KdCouple4     18
% 25: nCouple4
% 26: LatB const
% 27: Jasp mag const
% 28: Jasp sens const
% 29: KdLuc         19
% 30: MRTFRelease   20
% 31: KinsoloMRTF   21
% 32: Kin2MRTF      22
% 33: Kdim          
activeIdx = [1,3:7,9:13, 15,16,19,20,21,23,24,29,30,31,32];
paramNames = paramNames(activeIdx);
numParam = 22;
for i = 1:numParam
    InputOpts.Marginals(i).Type = 'Lognormal';
    InputOpts.Marginals(i).Name = paramNames{i};
    InputOpts.Marginals(i).Parameters = [0 log(4)];
end
myInput = uq_createInput(InputOpts);
SobolOpts.Type = 'uq_sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 100;
SobolOpts.Bootstrap.Replications = 10;
SobolOpts.Bootstrap.Alpha = 0.05;
SobolAnalysis = uq_createAnalysis(SobolOpts);
% MorrisSensOpts.Type = 'Sensitivity';
% MorrisSensOpts.Method = 'Morris';
% MorrisSensOpts.Morris.Cost = 1e4;
% MorrisAnalysis = uq_createAnalysis(MorrisSensOpts);

%% display
uq_visualizeSobolIndices(SobolAnalysis.Results, [1, 8, 15, 22])
% uq_display(MorrisAnalysis)