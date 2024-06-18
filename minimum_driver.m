%% Minimal test of mechano-circadian model
load('pSol_MAP.mat') % load pSol associated with maximum a posteriori parameters
stiffnessVals = [1,19,300,1e7]; % substrate stiffness values to test in kPa
timeSpan = [0, 3600*24*5]; % time span in seconds (5 days here)
figure
hold on
for i = 1:length(stiffnessVals)
    inhibVec = []; % empty is set to defaults, see MechanoCircadianModel for full description of inhibitor inputs
    [t,y, SSVals] = MechanoCircadianModel(timeSpan, stiffnessVals(i), pSol, inhibVec);
    plot(t/(3600*24), y(:,2))
end
xlabel('Time (days)')
ylabel('PER/CRY concentration (dimensionless)')
prettyGraph

%% Calculate period and amplitude as a function of substrate stiffness
stiffnessVals = logspace(-1, 3, 100); % log-spaced stiffness between 0.1 and 1000 kPa
period = zeros(size(stiffnessVals));
amplitude = zeros(size(stiffnessVals));
for i = 1:length(stiffnessVals)
    inhibVec = []; % empty is set to defaults, see MechanoCircadianModel for full description of inhibitor inputs
    [periodCur, amplitudeCur] = conditionToOutputs(pSol,stiffnessVals(i),inhibVec);
    period(i) = periodCur(2)/3600; % period for second variable (PER/CRY), converted to hours
    amplitude(i) = amplitudeCur(2); % dimensionless amplitude for PER/CRY
end
% plot period and relative amplitude in subplots
figure
subplot(2,1,1)
semilogx(stiffnessVals, period)
xlabel('Substrate stiffness (kPa)')
ylabel('Oscillation period (hr)')
prettyGraph
ylim([24 25.5])
subplot(2,1,2)
% normalize to amplitude on stiff substrate
semilogx(stiffnessVals, amplitude/amplitude(end))
xlabel('Substrate stiffness (kPa)')
ylabel('Normalized oscillation amplitude')
prettyGraph
ylim([0.5 2.5])