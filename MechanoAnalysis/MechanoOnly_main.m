%% Plot mechanotransduction dynamics with default parameters (for Fig 1B)
maxTime = 3600*6;
[T,Y] = MechanoOnlyModel([0 maxTime],[100,inf]);
figure
hold on
plot(T/60,Y(:,15)./(Y(:,17)+Y(:,18)),'LineWidth',1);
plot(T/60,Y(:,26)./Y(:,25),'LineWidth',1);
xlim([0 2*60])
legend('YAP/TAZ','MRTF')
ylabel('Nuclear to cytoplasmic ratio')
xlabel('Time (minutes)')
prettyGraph

%% extract pSolMRTF from full solve
p0 = [12*3600; 2; 1/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2];
if ~exist('myBayesianAnalysis','var')
    error('Load in myBayesianAnalysis first')
end
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP','burnIn',1000)
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
fixedParam = [2, 8, 14, 17, 22, 25];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
pSol = p0;
pSolMRTFIdx = [31, 32, 30, 18, 33, 27, 28, 19];
pSolMRTF = pSol(pSolMRTFIdx);
popParamMRTF = popParam(:, pSolMRTFIdx);

%% plot MRTF test cases
stiffnessVals = logspace(-1,4,100);
inhibVals = linspace(0,2,100);
MRTFEq = zeros(size(stiffnessVals));
YAPTAZEq = zeros(size(stiffnessVals));
FactinEq = zeros(size(stiffnessVals));
GactinEq = zeros(size(stiffnessVals));
for i = 1:length(stiffnessVals)
    inhibVec = [1,1,1,1, 0];
    inhibVec(6:9) = [1,1,3000,1];
    stiffnessVec = [stiffnessVals(i),inf];
    SSVar = MechanoSS(stiffnessVec, inhibVec, pSol);
    MRTFEq(i) = SSVar(25)/SSVar(26);
    YAPTAZEq(i) = SSVar(15)/(SSVar(17)+SSVar(18));
    FactinEq(i) = SSVar(5);
    GactinEq(i) = SSVar(9);
end
figure
semilogx(stiffnessVals,YAPTAZEq)
hold on
semilogx(stiffnessVals,MRTFEq)
hold on

%% plot cytD cases
cytDVals = 0:.1:12;
MRTFEq = zeros(size(cytDVals));
for i = 1:length(cytDVals)
    inhibVec = [1,1,1,1, cytDVals(i)];
    inhibVec(6:9) = [1,1,3000,1];
    stiffnessVec = [1e7,inf,0];
    SSVar = MechanoSS(stiffnessVec, inhibVec, pSol);
    MRTFEq(i) = SSVar(25)/SSVar(26);
    FActinEq(i) = SSVar(5);
    GActinEq(i) = SSVar(9);
    NPCAEq(i) = SSVar(8);
end
figure
plot(cytDVals,MRTFEq)
hold on

%% plot jasp cases
jaspVals = 0:.01:1.5;
MRTFEq = zeros(size(jaspVals));
for i = 1:length(jaspVals)
    actinPolyFactor = 1 + (1+pSolMRTF(5))*jaspVals(i)/pSolMRTF(6);
    inhibVec = [actinPolyFactor,1,1,1, 0];
    stiffnessVec = [1e7,inf,0];
    SSVar = MechanoSS(stiffnessVec, inhibVec, pSol);
    MRTFEq(i) = SSVar(25)/SSVar(26);
end
figure
plot(jaspVals,MRTFEq)
hold on

%% plot data comparison of MRTF vs. substrate stiffness (Fig S3)
stiffnessTests = [1e7; 1e7; 1e7; 0.175; 0.175; 0.175;...
        0.3; 20; 1e7;... 
        6.573705179; 15.79681275; 25.47144754; 34.58167331;...
        5; 100;...
        1.530002883; 13.77000253; 1.394647827; 14.82309684;...
        1e7; 1e7;...
        3; 40;...
        1e7; 1; 6; 50; 1e7; 1; 6; 50;...
        30; 30; 30; 0.3];
groups = {1:6, 7:9, 10:13, 14:15, 16:19, 20:21, 22:23, 24:31, 32:35};
MRTFEqData = [0.846153846; 6.692307692; 4.230769231; 0.615384615; 1.692307692;...
    0.923076923; 0.372093023; 2.511627907; 8; 1.651785714; 1.95610119; 1.943948413;...
    1.680272109; 1.923076923; 3.461538462; 1.201970443; 1.596059113; 1.448275862;...
    1.280788177; 5.68852459; 6.776818742; 0.64921466; 0.670157068; 2.798722045;...
    1.840255591; 1.878594249; 3.872204473; 2.485436893; 1.009708738; 1.048543689;...
    2.097087379; 1.345730564; 1.282941831; 1.254349864; 0.905891002]';
MRTFEqStd = [0.973008511; 3.447836128; 2.68239935; 0.486504255; 0.661717328;...
    0.407038663; 0.186046512; 0.465116279; 0.744186047; 0.144541983; 0.187825344; 0.130618877;...
    0.077764917; 0.307692308; 0.269230769; 0.313625496; 0.320919112; 0.393855274;...
    0.510553133; 0.163934426; 0.483353884; 0.217639892; 0.181366577; 0.93661217;...
    0.510879365; 0.397350617; 1.390727161; 0.919987005; 0.287495939; 0.402494315;...
    0.689990254; 0.078904541; 0.106530658; 0.075812277; 0.079320066]';
cytDIdx = [3,6,21,34];
jaspIdx = [2,5,33];
excludeIdx = [jaspIdx, cytDIdx, 22:23, 25:28];
keepLogic = true(length(MRTFEqStd),1);
keepLogic(excludeIdx) = false;
stiffnessVals = logspace(-1,4,100);
MRTFEq = zeros(size(popParam,1), length(stiffnessVals));
for k = 1:size(popParam,1)
    pCur = popParam(k,:);
    for i = 1:length(stiffnessVals)
        inhibVec = [1,1,1,1, 0];
        inhibVec(6:9) = [1,1,3000,1];
        stiffnessVec = [stiffnessVals(i),inf,0];
        SSVar = MechanoSS(stiffnessVec, inhibVec, pCur);
        MRTFEq(k,i) = SSVar(25)/SSVar(26);
    end
end

figure
colororder(linspecer(7))
markerList = {'o','*','x','s','d','','','p','^'};
for i = 1:length(groups)
    if i==7 || i==6
        continue
    end
    curIdx = groups{i};
    curIdx = curIdx(keepLogic(curIdx));
    errorbar(stiffnessTests(curIdx), MRTFEqData(curIdx), MRTFEqStd(curIdx),...
        MRTFEqStd(curIdx), 'LineStyle','none', 'LineWidth', 1, 'Marker', markerList{i})
    hold on
end
set(gca, 'XScale', 'log')
hold on
prctilePlot(stiffnessVals, MRTFEq')
xlim([.1 200])
ylim([0 4])
legend('McGee et al. 2011','Fearing et al. 2019','Hadden et al. 2017', 'Li et al. 2016',...
    'Hui et al. 2019', 'Sumey et al. 2023', 'Abenza et al. 2023','Model prediction',...
    '95% credible interval','Interquartile range')
prettyGraph
xlabel('Substrate stiffness (kPa)')
ylabel('MRTF nuclear to cytosolic ratio')

%% compile effects of cytoskeletal inhibitors (Fig S6)
cytDVals = 0:.1:10;
latBVals = 0:.05:5;
jaspVals = 0:.01:1;
cytDCell = {cytDVals, zeros(size(latBVals)), zeros(size(jaspVals))};
latBCell = {zeros(size(cytDVals)), latBVals, zeros(size(jaspVals))};
jaspCell = {zeros(size(cytDVals)), zeros(size(latBVals)), jaspVals};
conditionVecs = {cytDVals, latBVals, jaspVals};
treatmentStrings = {'Cytochalasin D (uM)','Latrunculin B (uM)','Jasplakinolide (uM)'};
figure
MRTFCell = cell(1,3);
YAPTAZCell = cell(1,3);
FActinCell = cell(1,3);
GActinCell = cell(1,3);
for i = 1:length(cytDCell)
    MRTFCell{i} = zeros(size(cytDCell{i}));
    YAPTAZCell{i} = zeros(size(cytDCell{i}));
    FActinCell{i} = zeros(size(cytDCell{i}));
    GActinCell{i} = zeros(size(cytDCell{i}));
    for j = 1:length(cytDCell{i})
        actinPolyFactor = 1 / (1 + (latBCell{i}(j)/pSol(26))) + (1 + pSol(27))*jaspCell{i}(j) / pSol(28);
        inhibVec = [actinPolyFactor,1,1,1, cytDCell{i}(j)];
        inhibVec(6:9) = [1,1,3000,1];
        stiffnessVec = [1e7,inf,0];
        SSVar = MechanoSS(stiffnessVec, inhibVec, pSol);
        MRTFCell{i}(j) = SSVar(25)/SSVar(26);
        YAPTAZCell{i}(j) = SSVar(15)/(SSVar(17)+SSVar(18));
        FActinCell{i}(j) = SSVar(5);
        GActinCell{i}(j) = SSVar(9);
    end
    subplot(3,2,2*i-1)
    curOrder = colororder("default");
    colororder(curOrder([1,4],:))
    plot(conditionVecs{i}, FActinCell{i}, 'LineWidth', 1.5, 'Color', curOrder(2,:))
    hold on
    plot(conditionVecs{i}, GActinCell{i}, 'LineWidth', 1.5, 'Color', curOrder(3,:))
    ylabel('Actin concentration')
    xlabel(treatmentStrings{i})
    if i==1
        legend('F-actin','G-actin')
    end
    ylim([0 510])
    prettyGraph
    subplot(3,2,2*i)
    plot(conditionVecs{i}, YAPTAZCell{i}, 'LineWidth', 1.5)
    ylabel('YAP/TAZ N/C')
    ylim([0 ceil(max(YAPTAZCell{i}))])
    yyaxis right
    plot(conditionVecs{i}, MRTFCell{i}, 'LineWidth', 1.5, 'Color', curOrder(4,:))
    ylabel('MRTF N/C')
    ylim([0 ceil(max(MRTFCell{i}))])
    xlabel(treatmentStrings{i})
    prettyGraph
    if i==1
        legend('YAP/TAZ N/C','MRTF N/C')
    end
end