bifFig = figure;
hold on
% omegaFig = figure;
% hold on
% omegaFig2 = figure;
% hold on
% stiffness = log(1e5);
% tauB1 = 8; pExpB = 2.5; KeB = 1;
% KiB = .04; KdBP = 0.5; KdB = 0.4;
% tauB2 = 8; pExpP = 2.5; KeP = 1;
% KaP = .5; KdP = .4;
% KeB2 = .1; KeP2 = .1;
% par = [tauB1*3600, pExpB, KeB/3600, KiB, KdBP/3600, KdB/3600,...
%     tauB2*3600, pExpP, KeP/3600, KaP, KdP/3600, .1/3600, .1/3600];
par = pSol;
par(29:30) = [2,0.1];
% par(12) = pSol(12) + pSol(23);
% par(13) = pSol(15) + pSol(20);
par(1) = par(1)/3600;
par(7) = par(7)/3600;
% par(20) = stiffness;
funcs=set_funcs(...
    'sys_rhs', @circadian_rhs,...
    'sys_tau',@()[1,7]);%tau is the first parameter in the parameter vector,...
%'sys_deri',@circadian_deri);
% KeP2Vals = pSol(13);%*[0,.1,.2,.5,1,2,5,10];
% KdBPVals = par(5)*[0,.1,.2,.5,1,2,5,10,50,100,500,1000];%[0,.1,.2,.5,1,2,5,10];%,5,10];
% KeB2Init = par(12)*[0,0,0,0,0,0,.5,1,2,2,2,2];%zeros(1,8);%[0,1,.1,.1,.1];
KdBPVals = par(5)*[0,.1,1,10];%[0,1,10,100,1000];
YInit = par(29)*[.5,.7,1,3];
% KeB2Init = par(12)*[0,0,.2,1.5];%[0,0,1,2,2];
plotLogic = false;
for i = 1:1%length(KdBPVals)
    par(5) = KdBPVals(i);
    par(29) = YInit(i);
%     par(13) = KeP2Vals(i);%.5;
    % compute ss
    stst.kind='stst';
    stst.parameter=par;
    stst.x = [0.1; 2];
    flag_newhheur=1; % flag_newhheur=1 is the default choice if this argument is omitted
    method=df_mthod(funcs,'stst',flag_newhheur);
    method.stability.minimal_real_part=-1;
    [stst,success]=p_correc(funcs,stst,[],[],method.point);
    % compute its stability:
    stst.stability=p_stabil(funcs,stst,method.stability);
    figure; clf;
    p_splot(stst); % plot its stability
    if ~plotLogic
        close(gcf)
    end

    figure; clf;
    % continue along branch
    ind_stiffness = 20;
%     ind_KeB2 = 12;
%     ind_KeP2 = 13;
    ind_Y = 29;
%     ind_M = 30;
    % get an empty branch with ind_KeB2 as a free parameter:
    branch1=df_brnch(funcs,ind_Y,'stst')
    branch1.parameter
    branch1.parameter.min_bound
    % set bounds for continuation parameter
    branch1.parameter.min_bound(1,:)=[ind_Y 0];
    branch1.parameter.max_bound(1,:)=[ind_Y 50];
    branch1.parameter.max_step(1,:)=[ind_Y .02];
    % use stst as a first branch point:
    branch1.point=stst;
    stst.parameter(ind_Y)=stst.parameter(ind_Y)+.02;
%     stst.parameter(ind_KeP2)=stst.parameter(ind_KeP2)+1e-5;
    [stst,success]=p_correc(funcs,stst,[],[],method.point)
    % use as a second branch point:
    branch1.point(2)=stst;
    % continue in one direction:
    [branch1,s,f,r]=br_contn(funcs,branch1,200)
%     % turn the branch around:
%     branch1=br_rvers(branch1);
%     % continue in the other direction:
%     [branch1,s,f,r]=br_contn(funcs,branch1,200)

    branch1.method.stability.minimal_real_part=-2;
    branch1=br_stabl(funcs,branch1,0,0);
    if ~plotLogic
        close(gcf)
    end

    % obtain suitable scalar measures to plot stability along branch:
    [xm,ym]=df_measr(1,branch1)
    figure; clf;
    br_plot(branch1,xm,ym,'b'); % plot stability along branch:
    ym.subfield='l0';
    br_plot(branch1,xm,ym,'c');
    % plot([0 5],[0 0],'-.');
    % axis([0 5 -2 1.5]);
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

    % Hopf bifurcation
    ind_hopf=find(arrayfun(@(x)real(x.stability.l0(1))>0,branch1.point),1,'last');
    hopf=p_tohopf(funcs,branch1.point(ind_hopf));
    method=df_mthod(funcs,'hopf',flag_newhheur); % get hopf calculation method parameters:
    method.stability.minimal_real_part=-1;
    [hopf,success]=p_correc(funcs,hopf,ind_Y,[],method.point) % correct hopf
    first_hopf=hopf;                    % store hopf point in other variable for later use
    hopf.stability=p_stabil(funcs,hopf,method.stability); % compute stability of hopf point
    figure; clf;
    p_splot(hopf);                     % plot stability of hopf point
    if ~plotLogic
        close(gcf)
    end

    % continue Hopf branch
    ind_Y = 29;
    ind_M = 30;
    branch2=df_brnch(funcs,[ind_Y,ind_M],'hopf'); % use hopf point as first point of hopf branch:
    branch2.parameter.min_bound(1:2,:)=[[ind_Y 0]' [ind_M 0]']';
    branch2.parameter.max_bound(1:2,:)=[[ind_Y 50]' [ind_M 50]']';
    branch2.parameter.max_step(1:2,:)=[[ind_Y .005]' [ind_M .005]']';
%     branch2=df_brnch(funcs,[ind_KeB2,ind_KdBP],'hopf'); % use hopf point as first point of hopf branch:
%     branch2.parameter.min_bound(1:2,:)=[[ind_KeB2 0]' [ind_KdBP 0]']';
%     branch2.parameter.max_bound(1:2,:)=[[ind_KeB2 1e-3]' [ind_KdBP 1e-3]']';
%     branch2.parameter.max_step(1:2,:)=[[ind_KeB2 1e-6]' [ind_KdBP 1e-6]']';
    branch2.point=hopf;

    hopf.parameter(ind_Y)=hopf.parameter(ind_Y)+5e-4; % perturb hopf point
%     hopf.parameter(ind_KdBP)=hopf.parameter(ind_KdBP)+.5e-5; % perturb hopf point
    [hopf,success]=p_correc(funcs,hopf,ind_Y,[],method.point); % correct hopf point, recompute stability
    branch2.point(2)=hopf;                                 % use as second point of hopf branch:
    figure; clf;
    [branch2,s,f,r]=br_contn(funcs,branch2,2000);            % continue with plotting hopf branch:
    branch2=br_rvers(branch2);                             % reverse Hopf branch
    [branch2,s,f,r]=br_contn(funcs,branch2,2000);            % continue in other direction
%     xlabel('KeB2');ylabel('KeP2');
    xlabel('Y');ylabel('M');
    if ~plotLogic
        close(gcf)
    end
    KdBPBranch = zeros(size(branch2.point));
    YBranch = zeros(size(branch2.point));
    MBranch = zeros(size(branch2.point));
    omegaBranch = zeros(size(branch2.point));
    for j = 1:length(branch2.point)
%         KdBPBranch(j) = branch2.point(j).parameter(5);
        YBranch(j) = branch2.point(j).parameter(29);
        MBranch(j) = branch2.point(j).parameter(30);
        omegaBranch(j) = branch2.point(j).omega;
    end
    figure(bifFig)
%     plot(KdBPBranch*3600,KeB2Branch*3600)
    plot(YBranch, MBranch)
%     figure(omegaFig)
%     plot(KeB2Branch*3600,2*pi./omegaBranch)
%     figure(omegaFig2)
%     plot(KeP2Branch*3600, 2*pi./omegaBranch)
end

%% continuation of periodic orbits
intervals=40;
degree=4;
[psol,stepcond]=p_topsol(funcs,first_hopf,1e-2,degree,intervals);
%% correct periodic solution guess:
method=df_mthod(funcs,'psol');
[psol,success]=p_correc(funcs,psol,4,stepcond,method.point)
% branch4=df_brnch(funcs,ind_stiffness,'psol'); % empty branch:
% branch4.parameter.min_bound(1,:)=[ind_stiffness -3];
% branch4.parameter.max_bound(1,:)=[ind_stiffness log(1e6)];
% branch4.parameter.max_step(1,:)=[ind_stiffness 1];
branch4=df_brnch(funcs,ind_KeB2,'psol'); % empty branch:
branch4.parameter.min_bound(1,:)=[ind_KeB2 0];
branch4.parameter.max_bound(1,:)=[ind_KeB2 .2e-3];
branch4.parameter.max_step(1,:)=[ind_KeB2 1e-6];
% make degenerate periodic solution with amplitude zero at hopf point:
deg_psol=p_topsol(funcs,first_hopf,0,degree,intervals);
% use deg_psol and psol as first two points on branch:
deg_psol.mesh=[];
branch4.point=deg_psol;
psol.mesh=[];
branch4.point(2)=psol;
figure; clf;
[branch4,s,f,r]=br_contn(funcs,branch4,50); % compute periodic solutions branch
xlabel('stiffness');ylabel('amplitude');
[xm,ym]=df_measr(0,branch4);
ym.field='period';
ym.col=1;
figure
br_plot(branch4,xm,ym,'b');% look at the period along the branch:

function y = circadian_rhs(xx,par)
    BLag1 = xx(1,2);
    BLag2 = xx(1,3);
    B = xx(1,1);
    P = xx(2,1);
%     SSVar = MechanoSS([exp(par(20)),inf]);
    % par = [tauB, nB, KeB, KiB, KdBP, KdB, tauP, nP, KeP, KaP, KdP,
    % KeB2max, K_YT, nY, KeP2max, K_MRTF, nM, cytoDConst, GActinThresh]
    nB = par(2);
    KeB = par(3)*3600;
    KiB = par(4);
    KdBP = par(5)*3600;
    KdB = par(6)*3600;
    nP = par(8);
    KeP = par(9)*3600;
    KaP = par(10);
    KdP = par(11)*3600;
    Y = par(29);
    M = par(30);
    KeB2 = 3600*(par(12)*Y^2/(par(13)^2+Y^2) + par(23)*M^2/(par(24)^2+M^2));
    KeP2 = 3600*(par(15)*M^2/(par(16)^2+M^2) + par(20)*Y^2/(par(21)^2+Y^2));
    y = [KeB/(1+(BLag1/KiB)^nB) + KeB2 - KdBP*B*P - KdB*B;
         KeP/(1+(KaP/BLag2)^nP) + KeP2 - KdBP*B*P - KdP*P];
end