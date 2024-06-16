%% Numerical bifurcation analysis for mechano-Circadian system
% this file was developed using code from the demo files for DDEBIFTOOL v3
% https://ddebiftool.sourceforge.net/demos/neuron/

bifFig = figure;
hold on
% load p0 from myBayesianAnalysis (Bayesian parameter est results)
p0 = [12*3600; 2; 0.01/3600; .04; 0.4/3600; 0.4/3600; 7.5*3600; 2; 0.1/3600; .5; 0.4/3600;...
    0.05/3600; 1; 2; .05/3600; 1; 2; 10; 3.25; .05/3600; 1; 2; .05/3600; 1; 2; 0.2; 2; 0.1; log(2)/(2*3600);...
    100; 1; 10; 2; 2; 1/3600; 0.1; 0.4/3600; 7.5*3600; .05/3600; 1; .05/3600; 1; 2; 2; 2];
if ~exist('myBayesianAnalysis','var')
    error('Load in myBayesianAnalysis first')
end
uq_postProcessInversionMCMC(myBayesianAnalysis,'PointEstimate','MAP','burnIn',500)
modeVals = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};
fixedParam = [5, 14, 17, 22, 25, 44, 45];
varyLogic = true(length(p0),1);
varyLogic(fixedParam) = false;
pSol = p0;
pSol(varyLogic) = pSol(varyLogic) .* modeVals';
par = zeros(1,length(pSol));
par(:) = pSol;
par(29:30) = [2, 2];
par(1) = par(1)/3600;
par(7) = par(7)/3600;
par(38) = par(38)/3600;
par([3,5,6,9,11,35,37]) = par([3,5,6,9,11,35,37])*3600; 
funcs=set_funcs(...
    'sys_rhs', @circadian_rhs,...
    'sys_tau',@()[1,7,38]);%tau is the first parameter in the parameter vector,...

nCoupleTest = false; % test different values of the Hill coefficient for mechanotransduction-circadian coupling
KdTest = true; % test different values for decay rates in circadian model
% define different parameter values to test
if KdTest
    KdBVals = par(6)*[1,0.1,3,10,1,1,1,1,1,1]; % KdB (p(6))
    KdPVals = par(11)*[1,1,1,1,.2,.5,2,1,1,1]; % KdP (p(11))
    KdRVals = par(37)*[1,1,1,1,1,1,1,.1,10,30]; % KdR (p(37))
    % initializations tailored to different parameter changes
    YInit = par(29)*[1.5,0.5,2,2,0.4,0.8,4,1.5,1.2,1.2];
    MInit = par(30)*[1.5,0.5,2,2,0.4,0.8,4,1.5,1.2,1.2];
    nMechVals = 2*ones(size(YInit)); %nMech
elseif nCoupleTest %#ok<UNRCH>
    nMechVals = [1,1.5,2,3,4.5];
    KdBVals = par(6)*ones(size(nMechVals)); % KdB (p(6))
    KdPVals = par(11)*ones(size(nMechVals)); % KdP (p(11))
    KdRVals = par(37)*ones(size(nMechVals)); % KdR (p(37))
    % initializations tailored to different parameter changes
    YInit = par(29)*1.5*ones(size(nMechVals));
    MInit = par(30)*1.5*ones(size(nMechVals));
end
plotLogic = false; % whether to leave all plots open
BranchesStored = cell(size(KdBVals)); % cell to store YAP/TAZ-MRTF bifurcation curves
for i = 1:length(KdBVals)
    par(6) = KdBVals(i);
    par(11) = KdPVals(i);
    par(37) = KdRVals(i);
    par(29) = YInit(i);
    par(30) = MInit(i);
    par(46) = nMechVals(i);
    
    % compute ss
    stst.kind='stst';
    stst.parameter=par;
    stst.x = [0.1; 2; 2];
    flag_newhheur=1; % flag_newhheur=1 is the default choice if this argument is omitted
    method=df_mthod(funcs,'stst',flag_newhheur);
    method.stability.minimal_real_part=-1;
    [stst,~]=p_correc(funcs,stst,[],[],method.point);
    % compute its stability:
    stst.stability=p_stabil(funcs,stst,method.stability);
    figure; clf;
    p_splot(stst); % plot its stability
    if ~plotLogic
        close(gcf)
    end

    figure; clf;
    % continue along branch with YAP/TAZ as a bifurcation parameter
    ind_Y = 29;
    branch1=df_brnch(funcs,ind_Y,'stst');
    % set bounds for continuation parameter
    branch1.parameter.min_bound(1,:)=[ind_Y 0];
    branch1.parameter.max_bound(1,:)=[ind_Y 50];
    branch1.parameter.max_step(1,:)=[ind_Y .02];
    % use stst as a first branch point:
    branch1.point=stst;
    stst.parameter(ind_Y)=stst.parameter(ind_Y)+.02; % perturb value of YAP/TAZ
    [stst,~]=p_correc(funcs,stst,[],[],method.point);
    % use as a second branch point:
    branch1.point(2)=stst;
    % continue in one direction:
    [branch1,~,~,~]=br_contn(funcs,branch1,200);
    branch1.method.stability.minimal_real_part=-2;
    branch1=br_stabl(funcs,branch1,0,0);
    if ~plotLogic
        close(gcf)
    end

    % obtain suitable scalar measures to plot stability along branch:
    [xm,ym]=df_measr(1,branch1);
    figure; clf;
    br_plot(branch1,xm,ym,'b'); % plot stability along branch:
    ym.subfield='l0';
    br_plot(branch1,xm,ym,'c');
    xlabel('Y');ylabel('\Re\lambda');
    if ~plotLogic
        close(gcf)
    end

    % plot stability versus point number:
    figure; clf;
    br_plot(branch1,[],ym,'b');
    br_plot(branch1,[],ym,'b.');
    plot([0 30],[0 0],'-.');
    xlabel('point number along branch');ylabel('\Re(\lambda)');
    if ~plotLogic
        close(gcf)
    end

    % Find the Hopf bifurcation point (only works if the above branch
    % crosses the Hopf bifurcation)
    ind_hopf=find(arrayfun(@(x)real(x.stability.l0(1))>0,branch1.point),1,'last');
    hopf=p_tohopf(funcs,branch1.point(ind_hopf));
    method=df_mthod(funcs,'hopf',flag_newhheur); % get hopf calculation method parameters:
    method.stability.minimal_real_part=-1;
    [hopf,~]=p_correc(funcs,hopf,ind_Y,[],method.point); % correct hopf
    first_hopf=hopf;                    % store hopf point in other variable for later use
    hopf.stability=p_stabil(funcs,hopf,method.stability); % compute stability of hopf point
    figure; clf;
    p_splot(hopf);                     % plot stability of hopf point
    if ~plotLogic
        close(gcf)
    end

    % continue Hopf branch by freeing up both YAP/TAZ and MRTF variables
    ind_Y = 29;
    ind_M = 30;
    branch2=df_brnch(funcs,[ind_Y,ind_M],'hopf'); % use hopf point as first point of hopf branch:
    branch2.parameter.min_bound(1:2,:)=[[ind_Y 0]' [ind_M 0]']';
    branch2.parameter.max_bound(1:2,:)=[[ind_Y 20]' [ind_M 20]']';
    branch2.parameter.max_step(1:2,:)=[[ind_Y .05]' [ind_M .05]']';
    branch2.point=hopf;
    store_hopf = hopf;

    store_hopf.parameter(ind_Y)=store_hopf.parameter(ind_Y)+5e-4; % perturb hopf point
    [hopf,~]=p_correc(funcs,store_hopf,ind_Y,[],method.point); % correct hopf point, recompute stability
    branch2.point(2)=hopf;                                 % use as second point of hopf branch:
    figure; clf;
    [branch2,~,~,~]=br_contn(funcs,branch2,2000);            % continue with plotting hopf branch:
    branch2=br_rvers(branch2);                             % reverse Hopf branch
    [branch2,~,~,~]=br_contn(funcs,branch2,2000);            % continue in other direction
    xlabel('Y');ylabel('M');
    if ~plotLogic
        close(gcf)
    end
    KdBPBranch = zeros(size(branch2.point));
    YBranch = zeros(size(branch2.point));
    MBranch = zeros(size(branch2.point));
    omegaBranch = zeros(size(branch2.point)); % save frequencies
    for j = 1:length(branch2.point)
        YBranch(j) = branch2.point(j).parameter(29);
        MBranch(j) = branch2.point(j).parameter(30);
        omegaBranch(j) = branch2.point(j).omega;
    end
    figure(bifFig)
    plot(YBranch, MBranch)
    BranchesStored{i} = {YBranch, MBranch};
end

% define right-hand-side (RHS) for circadian system (rates of change for
% each variable)
function y = circadian_rhs(xx,par)
    RLag1 = xx(3,2);
    BLag2 = xx(1,3);
    PLag2 = xx(2,3);
    BLag3 = xx(1,4);
    PLag3 = xx(2,4);
    B = xx(1,1);
    P = xx(2,1);
    R = xx(3,1);
    nB = par(2);
    KeB = par(3);
    KiB = par(4);
    % KdBP = par(5);
    KdB = par(6);
    nP0 = par(8);
    KeP = par(9);
    KaP = par(10);
    KdP = par(11);
    Y = par(29);
    M = par(30);
    KeP0_self = par(35);
    KiP = par(36);
    KdR = par(37);
    nP1 = par(43);
    CytoConv = 1.3851e6;
    NucConv = 3.3122e5;
    Ytot = 1.4784e6;
    Mtot = 1e6;
    YConc = Y*Ytot./(CytoConv + Y*NucConv);
    MConc = M*Mtot./(CytoConv + M*NucConv);
    nMech = par(46);
    KeB2 = 3600*(par(12)*YConc^nMech/(par(13)^nMech+YConc^nMech) +...
        par(23)*MConc^nMech/(par(24)^nMech+MConc^nMech));
    KeP2 = 3600*(par(15)*MConc^nMech/(par(16)^nMech+MConc^nMech) +...
        par(20)*YConc^nMech/(par(21)^nMech+YConc^nMech));
    KeR2 = 3600*(par(41)*MConc^nMech/(par(42)^nMech+MConc^nMech) +... 
        par(39)*YConc^nMech/(par(40)^nMech+YConc^nMech));
    y = [KeB/(1+(RLag1/KiB)^nB) + KeB2 - KdB*B;
         KeP/(1+(KaP/BLag2)^nP0) + KeP0_self/(1+(PLag2/KiP)^nP1) + KeP2 - KdP*P;
         KeP/(1+(KaP/BLag3)^nP0) + KeP0_self/(1+(PLag3/KiP)^nP1) + KeR2 - KdR*R];
end