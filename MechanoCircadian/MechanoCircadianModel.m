function [T,Y, SSVals] = MechanoCircadianModel(timeSpan, stiffnessParam, paramList, inhibVec, varargin)
% [T,Y,yinit,param] = EAFTests_Scott_Fraley_Rangamani_YAP_TAZ_model_2021__2D_and_3D_simulations_(argTimeSpan,argYinit,argParam)
%
% input:
%     argTimeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
%     argYinit is a vector of initial conditions for the state variables (optional)
%     argParam is a vector of values for the parameters (optional)
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yinit is the initial conditions that were used
%     param is the parameter vector that was used
%     allNames is the output solution variable names
%     allValues is the output solution variable values corresponding to the names
%
%     example of running this file: [T,Y,yinit,param,allNames,allValues] = myMatlabFunc; <-(your main function name)
%

if isempty(varargin)
    popVar = 0;
    noiseLevel = [0,0];
elseif length(varargin)==1
    popVar = varargin{1};
    noiseLevel = [0, 0];
elseif length(varargin)==2
    popVar = varargin{1};
    noiseLevel = varargin{2};
end
%
% Default Initial Conditions
%
yinit = [
	0.2;		% yinit(1) is the initial condition for 'CofilinP'
	0.7;		% yinit(2) is the initial condition for 'Fak'
	0.8;		% yinit(3) is the initial condition for 'mDia'
	0.0;		% yinit(4) is the initial condition for 'LaminA'
	17.9;		% yinit(5) is the initial condition for 'Fcyto'
	33.6;		% yinit(6) is the initial condition for 'RhoAGTP_MEM'
	0.0;		% yinit(7) is the initial condition for 'mDiaA'
	0.0;		% yinit(8) is the initial condition for 'NPCA'
	482.4;		% yinit(9) is the initial condition for 'Gactin'
	6.5;		% yinit(10) is the initial condition for 'NPC'
	0.0;		% yinit(11) is the initial condition for 'ROCKA'
	3.5;		% yinit(12) is the initial condition for 'Myo'
	1.8;		% yinit(13) is the initial condition for 'CofilinNP'
	3500.0;		% yinit(14) is the initial condition for 'LaminAp'
	0.7;		% yinit(15) is the initial condition for 'YAPTAZnuc'
	0.3;		% yinit(16) is the initial condition for 'Fakp'
	0.2;		% yinit(17) is the initial condition for 'YAPTAZP'
	0.7;		% yinit(18) is the initial condition for 'YAPTAZN'
	1.0;		% yinit(19) is the initial condition for 'RhoAGDP'
	1.9;		% yinit(20) is the initial condition for 'LIMK'
	1.5;		% yinit(21) is the initial condition for 'MyoA'
	1.0;		% yinit(22) is the initial condition for 'ROCK'
	1.0;		% yinit(23) is the initial condition for 'Positionboolean'
	0.1;		% yinit(24) is the initial condition for 'LIMKA'
];

%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%
param = [
	0.648;		% param(1) is 'krp'
	0.14;		% param(2) is 'kNC'
	0.4;		% param(3) is 'kra'
	0.7;		% param(4) is 'YAPTAZnuc_init_uM'
	0.015;		% param(5) is 'kf'
	4.0;		% param(6) is 'kmCof'
	4.0;		% param(7) is 'kfc1'
	0.625;		% param(8) is 'kdp'
	0.067;		% param(9) is 'kdmy'
	3.5;		% param(10) is 'Myo_init_uM'
	0.0;		% param(11) is 'mDiaA_init_uM'
	2.0;		% param(12) is 'kdl'
	1.5;		% param(13) is 'MyoA_init_uM'
	0.035;		% param(14) is 'kdf'
	8.7;		% param(15) is 'krNPC'
	0.0;		% param(16) is 'ROCKA_init_uM'
	25.0;		% param(17) is 'Emol'
	0.34;		% param(18) is 'kcatcof'
	1.8255;		% param(19) is 'SAV'
	0.8;		% param(20) is 'kdrock'
	1.0;		% param(21) is 'RhoAGDP_init_uM'
	0.0168;		% param(22) is 'kfkp'
	2300.0;		% param(23) is 'Size_Cyto'
	0.005;		% param(24) is 'kdmdia'
	50.0;		% param(25) is 'alpha'
	1704.7954979572205;		% param(26) is 'Size_ECM'
	0.0;		% param(27) is 'Voltage_PM'
	0.03;		% param(28) is 'kmr'
	77.56;		% param(29) is 'gamma'
	0.002;		% param(30) is 'kmp'
	1.8;		% param(31) is 'CofilinNP_init_uM'
	16.0;		% param(32) is 'k11'
	1.0;		% param(33) is 'netValence_r5'
	1.0;		% param(34) is 'netValence_r4'
	1.0;		% param(35) is 'netValence_r3'
	1.0;		% param(36) is 'netValence_r1'
	0.46;		% param(37) is 'kflaminA'
	120.0;		% param(38) is 'kly'
	0.3;		% param(39) is 'Fakp_init_uM'
	0.07;		% param(40) is 'klr'
	16.0;		% param(41) is 'kll'
	0.0;		% param(42) is 'LaminA_init_molecules_um_2'
	0.165;		% param(43) is 'mDiaB'
	0.0;		% param(44) is 'Voltage_NM'
	33.6;		% param(45) is 'RhoAGTP_MEM_init_molecules_um_2'
	2.8E-7;		% param(46) is 'kfNPC'
	9.64853321E-5;		% param(47) is 'mlabfix_F_nmol_'
	602.2;		% param(48) is 'unitconversionfactor'
	300.0;		% param(49) is 'mlabfix_T_'
	0.8;		% param(50) is 'mDia_init_uM'
	1000.0;		% param(51) is 'K_millivolts_per_volt'
	0.001;		% param(52) is 'Kr_r15'
	0.0;		% param(53) is 'Kr_r12'
	100.0;		% param(54) is 'Clamin'
	36.0;		% param(55) is 'epsilon'
	0.2;		% param(56) is 'YAPTAZP_init_uM'
	17.9;		% param(57) is 'Fcyto_init_uM'
	0.3;		% param(58) is 'ROCKB'
	482.4;		% param(59) is 'Gactin_init_uM'
	1260.0;		% param(60) is 'Size_PM'
	0.1;		% param(61) is 'LIMKA_init_uM'
	550.0;		% param(62) is 'Size_Nuc'
	10.0;		% param(63) is 'kin2'
	3.141592653589793;		% param(64) is 'mlabfix_PI_'
	9.0E-6;		% param(65) is 'p'
	96485.3321;		% param(66) is 'mlabfix_F_'
	1.0;		% param(67) is 'kinSolo2'
	8314.46261815;		% param(68) is 'mlabfix_R_'
	55.49;		% param(69) is 'tau'
	1.0;		% param(70) is 'kout2'
	1.0;		% param(71) is 'ROCK_init_uM'
	1.0E-9;		% param(72) is 'mlabfix_K_GHK_'
	3.5;		% param(73) is 'kdep'
	393.0;		% param(74) is 'Size_NM'
	0.0;		% param(75) is 'Kr_r7'
	0.0;		% param(76) is 'Kr_r6'
	0.0;		% param(77) is 'Kr_r5'
	3.25;		% param(78) is 'C'
	0.0;		% param(79) is 'Kr_r4'
	0.0;		% param(80) is 'Kr_r3'
	0.0;		% param(81) is 'Kr_r2'
	0.0;		% param(82) is 'Kr_r1'
	0.0;		% param(83) is 'Kr_r0'
	6.5;		% param(84) is 'NPC_init_molecules_um_2'
	6.02214179E11;		% param(85) is 'mlabfix_N_pmol_'
	7.6E-4;		% param(86) is 'kCY'
	0.2;		% param(87) is 'CofilinP_init_uM'
	3500.0;		% param(88) is 'LaminAp_init_molecules_um_2'
	0.56;		% param(89) is 'kCN'
	1.0;		% param(90) is 'netValence_r16'
	1.0;		% param(91) is 'netValence_r15'
	1.0;		% param(92) is 'netValence_r14'
	1.9;		% param(93) is 'LIMK_init_uM'
	1.0;		% param(94) is 'netValence_r12'
	0.04;		% param(95) is 'kturnover'
	0.7;		% param(96) is 'Fak_init_uM'
	0.379;		% param(97) is 'ksf'
	0.7;		% param(98) is 'YAPTAZN_init_uM'
	5.0;		% param(99) is 'n2'
	2.6;		% param(100) is 'n1'
	0.0;		% param(101) is 'NPCA_init_molecules_um_2'
	1.0;		% param(102) is 'Positionboolean_init_molecules_um_2'
	0.001660538783162726;		% param(103) is 'KMOLE'
];
param(17) = stiffnessParam(1); % substrate stiffness
param(104) = stiffnessParam(2); % time parameter for stiffness increase

if length(popVar) == 108
    varVecCur = popVar(1:104);
    varVecCur = varVecCur(:);
else
    varVecCur = sqrt(popVar).*randn(size(param(1:104)));
end
varVecCur([33:36,47:49,51,64,66,68,72,85,90:92,94,99:100,102:103]) = 0;
param(1:104) = param(1:104) .* exp(varVecCur);


param(105:121) = paramList(1:17);
% load inhibition parameters [actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)]
param(3) = param(3)*inhibVec(1); % inhibition of actin polymerization
param([55,69]) = param([55,69])*inhibVec(2); % inhibit ROCK-mediated catalysis
param(119) = param(119)*inhibVec(3); % inhibit coupling between MRTF and circadian
foldYAPOverexpress = inhibVec(4);
% param(2) = param(2)*inhibVec(4); % inhibit YAP phosphorylation
% if inhibVec(4) > 0
    % param([4,56,98]) = 2*param([4,56,98]);
% end
% if applicable, assess cytoD conc
param(123) = 0;
if length(inhibVec)>4
    param(123) = inhibVec(5);
end
if length(inhibVec)>5
    param(2) = inhibVec(6)*param(2);
end
if length(inhibVec)>6
    param(86) = inhibVec(7)*param(86); % inhibit actin-myosin interactions (blebbistatin) - YAPTAZ dephos
    param(46) = inhibVec(7)*param(46); % inhibit actin-myosin interactions (blebbistatin) - NPC opening
    param([5,97]) = inhibVec(8)*param([5,97])/3000; % account for contact area
    param(52) = inhibVec(9)*param(52); % inhibit phosph. of LaminA
end
param(78) = paramList(19);


param(122) = paramList(30);%172.93; %43.2495; %112.5; %MRTFReleaseConst
param(124) = 1e6;%3.15e6; %1.0724e7; %1.0153; %MRTFTot
param(125) = paramList(31);%6.289; %8.4028; %2.6383; %Kinsolo_MRTF
param(126) = paramList(32);%53.98; %16.5570; %6.2886; %Kin2_MRTF
% param(127) = 0.8232; %0.1346; %0.2725; %kout_MRTF
param(127) = paramList(34);%0.461; %0.3383; %0.25; %Kcap (cytoD)
param(128) = paramList(33); %Kdim (cytoD)
MRTFVariableParam = [122,124,125,126];
if length(popVar) == 108
    varVecMRTF = popVar(105:108);
else
    varVecMRTF = sqrt(popVar).*randn(size(MRTFVariableParam));
end
varVecMRTF = varVecMRTF(:);
param(MRTFVariableParam) = param(MRTFVariableParam) .* exp(varVecMRTF);

% YAP/TAZ-PER/CRY coupling
param(129:131) = paramList(20:22);
% MRTF-BMAL1 coupling
param(132:134) = paramList(23:25);
param(135) = paramList(29);

% 5SA-YAP overexpression
param(136) = foldYAPOverexpress;
%
% invoke the integrator
%
% [T,Y] = ode15s(@f,timeSpan,yinit,odeset('RelTol',1e-8),param,yinit);
% [TRed,YRed] = ode15s(@fRed, timeSpan, yinit([4,15]), odeset('RelTol',1e-8),param);
% [SSVals,tauVals] = fSS(param);
% SSVarMat = zeros(length(TRed),24);
% for k = 1:length(TRed)
%     [~,SSVarMat(k,:)] = fRed(TRed(k),YRed(k,:),param);
% end
[T,Y] = fDDE(timeSpan,param,noiseLevel);
SSVals = fSS(param);

end


% reduced dde
function [t,y] = fDDE(timeSpan, p, noiseLevel)
	% State Variables

	% Constants
	krp = p(1);
	kNC = p(2);
	kra = p(3);
	YAPTAZnuc_init_uM = p(4);
	kf = p(5);
	kmCof = p(6);
	kfc1 = p(7);
	kdp = p(8);
	kdmy = p(9);
	Myo_init_uM = p(10);
	mDiaA_init_uM = p(11);
	kdl = p(12);
	MyoA_init_uM = p(13);
	kdf = p(14);
	krNPC = p(15);
	ROCKA_init_uM = p(16);
	kcatcof = p(18);
	SAV = p(19);
	kdrock = p(20);
	RhoAGDP_init_uM = p(21);
	kfkp = p(22);
	Size_Cyto = p(23);
	kdmdia = p(24);
	alpha = p(25);
	Size_ECM = p(26);
	Voltage_PM = p(27);
	kmr = p(28);
	gamma = p(29);
	kmp = p(30);
	CofilinNP_init_uM = p(31);
	k11 = p(32);
	netValence_r5 = p(33);
	netValence_r4 = p(34);
	netValence_r3 = p(35);
	netValence_r1 = p(36);
	kflaminA = p(37);
	kly = p(38);
	Fakp_init_uM = p(39);
	klr = p(40);
	kll = p(41);
	LaminA_init_molecules_um_2 = p(42);
	mDiaB = p(43);
	Voltage_NM = p(44);
	RhoAGTP_MEM_init_molecules_um_2 = p(45);
	kfNPC = p(46);
	mlabfix_F_nmol_ = p(47);
	unitconversionfactor = p(48);
	mlabfix_T_ = p(49);
	mDia_init_uM = p(50);
	K_millivolts_per_volt = p(51);
	Kr_r15 = p(52);
	Kr_r12 = p(53);
	Clamin = p(54);
	epsilon = p(55);
	YAPTAZP_init_uM = p(56);
	Fcyto_init_uM = p(57);
	ROCKB = p(58);
	Gactin_init_uM = p(59);
	Size_PM = p(60);
	LIMKA_init_uM = p(61);
	Size_Nuc = p(62);
	kin2 = p(63);
	mlabfix_PI_ = p(64);
	propStiff = p(65);
	mlabfix_F_ = p(66);
	kinSolo2 = p(67);
	mlabfix_R_ = p(68);
	tau = p(69);
	kout2 = p(70);
	ROCK_init_uM = p(71);
	mlabfix_K_GHK_ = p(72);
	kdep = p(73);
	Size_NM = p(74);
	Kr_r7 = p(75);
	Kr_r6 = p(76);
	Kr_r5 = p(77);
	C = p(78);
	Kr_r4 = p(79);
	Kr_r3 = p(80);
	Kr_r2 = p(81);
	Kr_r1 = p(82);
	Kr_r0 = p(83);
	NPC_init_molecules_um_2 = p(84);
	mlabfix_N_pmol_ = p(85);
	kCY = p(86);
	CofilinP_init_uM = p(87);
	LaminAp_init_molecules_um_2 = p(88);
	kCN = p(89);
	netValence_r16 = p(90);
	netValence_r15 = p(91);
	netValence_r14 = p(92);
	LIMK_init_uM = p(93);
	netValence_r12 = p(94);
	kturnover = p(95);
	Fak_init_uM = p(96);
	ksf = p(97);
	YAPTAZN_init_uM = p(98);
	n2 = p(99);
	n1 = p(100);
	NPCA_init_molecules_um_2 = p(101);
	Positionboolean_init_molecules_um_2 = p(102);
	KMOLE = p(103);
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);
    
    %Circadian parameters
    tauB = p(105);
    pExpB = p(106);
    KeB = p(107);
    KiB = p(108);
    KdBP = p(109);
    KdB = p(110);
    tauP = p(111);
    pExpP = p(112);
    KeP = p(113);
    KaP = p(114);
    KdP = p(115);

    magCouple1 = p(116);
    KdCouple1 = p(117);
    nCouple1 = p(118);
    magCouple2 = p(119);
    KdCouple2 = p(120);
    nCouple2 = p(121);

    MRTFReleaseConst = p(122);
    cytoDConc = p(123);

    MRTFTot = p(124);
    Kinsolo_MRTF = p(125);
    Kin2_MRTF  = p(126);
    Kcap = p(127);
    Kdim = p(128);

    magCouple3 = p(129); % for YAP/TAZ coupling to PER/CRY
    KdCouple3 = p(130);
    nCouple3 = p(131);
    magCouple4 = p(132); % for MRTF coupling to BMAL1
    KdCouple4 = p(133);
    nCouple4 = p(134);

    KdLuc = p(135);

    
    lags = [tauB, tauP, 12*3600]; % third is nominal delay (only affects time shift of reporter, can be any delay in general)

%     DDESol = dde23(@ddefun, lags, @history, timeSpan);
    SSVals = fSS(p);
    YAPTAZnuc_SS = SSVals(15);
    MRTFnuc_SS = SSVals(25);
    KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
           magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4);
    KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
           magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3);
    fsolveOptions = optimoptions('fsolve','Display','off');
    unstableSS = fsolve(@(ss) [KeB/(1+(ss(1)/KiB)^pExpB) + KeB2 - KdBP*ss(1)*ss(2) - KdB*ss(1);
                               KeP/(1+(KaP/ss(1))^pExpP) + KeP2 - KdBP*ss(1)*ss(2) - KdP*ss(2)],...
                [0.1,0.1], fsolveOptions);
    unstableSS(3) = 0;
    
    startTime = tic;
    if all(noiseLevel == 0)
        try DDESol = dde23(@ddefun_YAPSS, lags, @history_CircOnly, timeSpan);
            t = DDESol.x';
            y = DDESol.y';
            [~, lucInt] = ode15s(@odeLuc, t, 0, odeset(), t, y(:,1));
            y(:,3) = lucInt';
        catch
            t = (timeSpan(1):3600:timeSpan(end))';
            y = zeros(length(t), 3);
        end
    elseif all(noiseLevel > 0)
        [t, y] = sdde_circadian_solve(noiseLevel(1), noiseLevel(2), noiseLevel(3), lags);
    else
        error('Noise level must be greater than or equal to zero')
    end
	
    function dydt = ddefun(t,y,Z)
        clockLag1 = Z(4,1);
        clockLag2 = Z(4,2);
        LaminA = y(1);
        YAPTAZnuc = y(2);
        MRTFnuc = y(3);

%         if t < 48*3600
%             YAPTAZTime = 0;
%         else
%             YAPTAZTime = t - 48*3600;
%         end
        YAPTAZTime = t;
        if isinf(p(104))
            Emol = p(17);
        else
            Emol = p(17)*YAPTAZTime/p(104);
        end

        %SS calcs (some analytical)
        Faktot = Fak_init_uM + Fakp_init_uM;
        Fakp = Faktot * (ksf*Emol + kf *(C+Emol)) / ((kf+kdf)*(C + Emol) + ksf*Emol);
        Fak = Faktot - Fakp;

        RhoAGTP_init_uM = RhoAGTP_MEM_init_molecules_um_2 * (UnitFactor_uM_um3_molecules_neg_1*KFlux_PM_Cyto);
        RhoATot = RhoAGDP_init_uM + RhoAGTP_init_uM;
        RhoAGTP = RhoATot * kfkp*(gamma*Fakp^n2 + 1) / (kdp + kfkp*(gamma*Fakp^n2 + 1));
        RhoAGDP = RhoATot - RhoAGTP;
        RhoAGTP_MEM = RhoAGTP * (unitconversionfactor * SAV);

        ROCKTot = ROCKA_init_uM + ROCK_init_uM;
        ROCKAInf = ROCKTot * (krp*RhoAGTP/(kdrock + krp*RhoAGTP));
        ROCKA0 = ROCKA_init_uM;
        ROCKA = ROCKAInf - (ROCKAInf-ROCKA0)*exp(-(kdrock + krp*RhoAGTP)*YAPTAZTime);
%         ROCKA = ROCKAInf;
        ROCK = ROCKTot - ROCKA;

        mDiaTot = mDia_init_uM + mDiaA_init_uM;
        mInf = mDiaTot*kmp*RhoAGTP/(kmp*RhoAGTP + kdmdia);
        m0 = mDiaA_init_uM;
        mDiaA = mInf - (mInf-m0)*exp(-(kmp*RhoAGTP + kdmdia)*YAPTAZTime);
        %     mDiaA = mInf;
        mDia = mDiaTot - mDiaA;

        % Compute other SS values
        % define smooth/tanh functions
        smoothROCKA = (ROCKA/2) * (tanh(20*(ROCKA-ROCKB)) + 1);
        smoothmDiaA = (mDiaA/2) * (tanh(20*(mDiaA-mDiaB)) + 1);
        LIMKTot = LIMKA_init_uM + LIMK_init_uM;
        LIMKA = LIMKTot*klr*(tau*smoothROCKA + 1)/(kdl + klr*(tau*smoothROCKA + 1));
        LIMK = LIMKTot - LIMKA;
        CofilinTot = CofilinP_init_uM + CofilinNP_init_uM;
        CofilinRoots = roots([kturnover, ...
            kcatcof*LIMKA + kmCof*kturnover - kturnover*CofilinTot, ...
            -kturnover*kmCof*CofilinTot]);
        CofilinNPinf = CofilinRoots(CofilinRoots > 0);
        CofilinNP = CofilinNPinf;
        %     CofilinB = (kturnover + kcatcof*LIMKA/kmCof);
        %     CofilinNP = CofilinNPinf - (CofilinNPinf-CofilinNP_init_uM)*exp(-CofilinB*t);
        CofilinP = CofilinTot - CofilinNP;

        ActinTot = Fcyto_init_uM + Gactin_init_uM;
        Fcyto = ActinTot*kra*(alpha*smoothmDiaA+1)...
            /(kdep + kfc1*CofilinNP + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/cytoDConst));
        Gactin = ActinTot - Fcyto*(1+cytoDConc/cytoDConst);

        MyoTot = Myo_init_uM + MyoA_init_uM;
        MyoA = MyoTot*kmr*(epsilon*smoothROCKA+1)/(kdmy + kmr*(epsilon*smoothROCKA+1));
        Myo = MyoTot - MyoA;

        % NPC from SS
        Ecytosol = propStiff*(Fcyto^n1);
        LaminATot = LaminAp_init_molecules_um_2 + LaminA_init_molecules_um_2;
        %     LaminA = LaminATot*kflaminA*Ecytosol/(kflaminA*Ecytosol + Kr_r15*(Clamin+Ecytosol));
        LaminAp = LaminATot - LaminA;
        NPCTot = NPCA_init_molecules_um_2 + NPC_init_molecules_um_2;
        NPCA = NPCTot*kfNPC*LaminA*Fcyto*MyoA / (krNPC + kfNPC*LaminA*Fcyto*MyoA);
        NPC = NPCTot - NPCA;

        % YAPTAZ (calc in terms of molecules, then convert)
        CytoConvert = Size_Cyto/UnitFactor_uM_um3_molecules_neg_1;
        NucConvert = Size_Nuc/UnitFactor_uM_um3_molecules_neg_1;
        YAPTAZnuc = YAPTAZnuc*NucConvert; % convert to molecules
        YAPTAZTot =  (YAPTAZP_init_uM*CytoConvert + YAPTAZN_init_uM*CytoConvert ...
            + YAPTAZnuc_init_uM*NucConvert);
        YAPTAZPFraction = (kNC/(kNC + kCN + kCY*Fcyto*MyoA)); % rapid phosphorylation
        %     YAPTAZnucFraction = (kinSolo2 + kin2*NPCA)/(kout2*Size_Cyto/Size_Nuc + kinSolo2 + kin2*NPCA);
        YAPTAZCytoTot = YAPTAZTot - YAPTAZnuc;
        YAPTAZP = YAPTAZCytoTot * YAPTAZPFraction;
        YAPTAZN = YAPTAZCytoTot - YAPTAZP;
        YAPTAZNPTot = YAPTAZTot - YAPTAZP;
        %     YAPTAZCytoTot = YAPTAZTot / (1 + (1-YAPTAZPFraction)*YAPTAZnucFraction/(1-YAPTAZnucFraction));
        %     YAPTAZP = YAPTAZCytoTot * YAPTAZPFraction;
        %     YAPTAZN = YAPTAZCytoTot - YAPTAZP;
        %     YAPTAZnuc = YAPTAZTot - YAPTAZCytoTot;
        YAPTAZN = YAPTAZN/CytoConvert;
        YAPTAZP = YAPTAZP/CytoConvert;
        YAPTAZnuc = YAPTAZnuc/NucConvert;

        %MRTF
        MRTFFactor = 1/(1+(Gactin/MRTFReleaseConst)^2);
        nucPerm_MRTF = kinsolo_MRTF + kin2_MRTF*NPCA;

        nucPerm = kinSolo2+kin2*NPCTot*kfNPC*y(1)*Fcyto*MyoA/(kfNPC*y(1)*Fcyto*MyoA + krNPC);
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3);

        % Rates (LaminA, YAPTAZ, MRTF, BMAL, and PER dynamics)
        if YAPTAZTime <= 0
            dydt = [0; 0; 0; 
                    KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdBP*y(4)*y(5) - KdB*y(4) + KeB2;
                    KeP*(1./(1+(KaP/clockLag2).^pExpP)) - KdBP*y(4)*y(5) - KdP*y(5) + KeP2];
        else
            dydt = [kflaminA*(Ecytosol/(Ecytosol+Clamin))*(LaminATot - y(1)) - Kr_r15*y(1);
                Size_NM*(nucPerm*((YAPTAZTot-y(2)*NucConvert)*(1-YAPTAZPFraction)/CytoConvert) - kout2*y(2))/NucConvert;
                Size_NM*(nucPerm_MRTF*((MRTFTot-y(3)*NucConvert)*MRTFFactor/CytoConvert) - kout_MRTF*y(3))/NucConvert;
                KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdBP*y(4)*y(5) - KdB*y(4) + KeB2;
                KeP*(1./(1+(KaP/clockLag2).^pExpP)) - KdBP*y(4)*y(5) - KdP*y(5) + KeP2];
        end
    end

    function s = history(t)
        s = [LaminA_init_molecules_um_2;
             YAPTAZnuc_init_uM;
             YAPTAZnuc_init_uM;
             0;
             0];
    end

    function dydt = ddefun_YAPSS(t,y,Z)
        clockLag1 = Z(1,1);
        clockLag2 = Z(1,2);
        mechanoStartEffect = false;
        if t < lags(1) && mechanoStartEffect
            YAPTAZnuc_B = YAPTAZnuc_init_uM;
            MRTFnuc_B = .301;
        else
            YAPTAZnuc_B = YAPTAZnuc_SS;
            MRTFnuc_B = MRTFnuc_SS;
        end
        if t < lags(2) && mechanoStartEffect
            YAPTAZnuc_P = YAPTAZnuc_init_uM;
            MRTFnuc_P = .301;
        else
            YAPTAZnuc_P = YAPTAZnuc_SS;
            MRTFnuc_P = MRTFnuc_SS;
        end
        % clockLag3 = Z(1,3);
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_B)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_B)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_P)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_P)^nCouple3);

        % Rates (BMAL and PER dynamics)
        dydt = [KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdBP*y(1)*y(2) - KdB*y(1) + KeB2;
                KeP*(1./(1+(KaP/clockLag2).^pExpP)) - KdBP*y(1)*y(2) - KdP*y(2) + KeP2];
                % KeP*(1./(1+(KaP/clockLag3).^pExpP)) - KdLuc*y(3) + KeP2];
        curEval = toc(startTime);
        if curEval > 10
            error('Solver took too long')
        end
    end

    function s = history_CircOnly(t)
        s = [5*unstableSS(1); 0.2*unstableSS(2)];%; unstableSS(3)];
    end

    function dluc = odeLuc(t, luc, tStored, BStored)
        BCur = interp1(tStored, BStored, t);
        dluc = KeP*(1./(1+(KaP/BCur).^pExpP)) - KdLuc*luc + KeP2;
    end

    function [t,y] = sdde_circadian_solve(eps_B, eps_P, eps_L, lags)

        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1) +...
               magCouple4 / (1 + (KdCouple4/MRTFnuc_SS)^nCouple4);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2) +...
               magCouple3 / (1 + (KdCouple3/YAPTAZnuc_SS)^nCouple3);
        dt = 2; % ~2 s time step
        delayFactor1 = round(lags(1)/dt);
        delayFactor2 = round(lags(2)/dt);
        t = (0:dt:timeSpan(2))';
        y = zeros(length(t),3);
        y(1,1:2) = history_CircOnly(0);
        y(1,3) = 0;
        dy = zeros(1,3);
        eta_B = 0;
        eta_P = 0;
        eta_L = 0;
        rand_stored = sqrt(dt)*randn(length(t)-1, 3);
        for i = 2:length(t)
            eta_B = eta_B - dt*eta_B + eps_B*rand_stored(i-1,1) + 0.5*eps_B^2*(rand_stored(i-1,1)^2-dt);
            eta_P = eta_P - dt*eta_P + eps_P*rand_stored(i-1,2) + 0.5*eps_P^2*(rand_stored(i-1,2)^2-dt);
            eta_L = eta_L - dt*eta_L + eps_L*rand_stored(i-1,3) + 0.5*eps_L^2*(rand_stored(i-1,3)^2-dt);
            % Rates (BMAL and PER dynamics)
            B = y(i-1,1);
            P = y(i-1,2);
            L = y(i-1,3);
            if t(i-1) < lags(1)
                BLag1 = y(1,1);
            else
                BLag1 = y(i-1-delayFactor1,1);
            end
            if t(i-1) < lags(2)
                BLag2 = y(1,1);
            else
                BLag2 = y(i-1-delayFactor2,1);
            end
            dy(1) = KeB*(1./(1+(BLag1/KiB).^pExpB)) - KdBP*B*P - KdB*B + KeB2 + eta_B;
            dy(2) = KeP*(BLag2^pExpP./(BLag2^pExpP+KaP^pExpP)) - KdBP*B*P - KdP*P + KeP2 + eta_P;
            dy(3) = KeP*(1./(1+(KaP/B).^pExpP)) - KdLuc*L + KeP2 + eta_L;
            y(i,:) = y(i-1,:) + dy*dt;
            if any(y(i,:) < 0)
                yCur = y(i,:);
                yCur(yCur < 0) = 0;
                y(i,:) = yCur;
            end
        end
    end
end


function [SSVar, tauVals] = fSS(p)
	% Constants
	krp = p(1);
	kNC = p(2);
	kra = p(3);
	YAPTAZnuc_init_uM = p(4);
	kf = p(5);
	kmCof = p(6);
	kfc1 = p(7);
	kdp = p(8);
	kdmy = p(9);
	Myo_init_uM = p(10);
	mDiaA_init_uM = p(11);
	kdl = p(12);
	MyoA_init_uM = p(13);
	kdf = p(14);
	krNPC = p(15);
	ROCKA_init_uM = p(16);
	Emol = p(17);
%     Emol = Emol*(t/(3600*48));
	kcatcof = p(18);
	SAV = p(19);
	kdrock = p(20);
	RhoAGDP_init_uM = p(21);
	kfkp = p(22);
	Size_Cyto = p(23);
	kdmdia = p(24);
	alpha = p(25);
	Size_ECM = p(26);
	Voltage_PM = p(27);
	kmr = p(28);
	gamma = p(29);
	kmp = p(30);
	CofilinNP_init_uM = p(31);
	k11 = p(32);
	netValence_r5 = p(33);
	netValence_r4 = p(34);
	netValence_r3 = p(35);
	netValence_r1 = p(36);
	kflaminA = p(37);
	kly = p(38);
	Fakp_init_uM = p(39);
	klr = p(40);
	kll = p(41);
	LaminA_init_molecules_um_2 = p(42);
	mDiaB = p(43);
	Voltage_NM = p(44);
	RhoAGTP_MEM_init_molecules_um_2 = p(45);
	kfNPC = p(46);
	mlabfix_F_nmol_ = p(47);
	unitconversionfactor = p(48);
	mlabfix_T_ = p(49);
	mDia_init_uM = p(50);
	K_millivolts_per_volt = p(51);
	Kr_r15 = p(52);
	Kr_r12 = p(53);
	Clamin = p(54);
	epsilon = p(55);
	YAPTAZP_init_uM = p(56);
	Fcyto_init_uM = p(57);
	ROCKB = p(58);
	Gactin_init_uM = p(59);
	Size_PM = p(60);
	LIMKA_init_uM = p(61);
	Size_Nuc = p(62);
	kin2 = p(63);
	mlabfix_PI_ = p(64);
	propStiff = p(65);
	mlabfix_F_ = p(66);
	kinSolo2 = p(67);
	mlabfix_R_ = p(68);
	tau = p(69);
	kout2 = p(70);
	ROCK_init_uM = p(71);
	mlabfix_K_GHK_ = p(72);
	kdep = p(73);
	Size_NM = p(74);
	Kr_r7 = p(75);
	Kr_r6 = p(76);
	Kr_r5 = p(77);
	C = p(78);
	Kr_r4 = p(79);
	Kr_r3 = p(80);
	Kr_r2 = p(81);
	Kr_r1 = p(82);
	Kr_r0 = p(83);
	NPC_init_molecules_um_2 = p(84);
	mlabfix_N_pmol_ = p(85);
	kCY = p(86);
	CofilinP_init_uM = p(87);
	LaminAp_init_molecules_um_2 = p(88);
	kCN = p(89);
	netValence_r16 = p(90);
	netValence_r15 = p(91);
	netValence_r14 = p(92);
	LIMK_init_uM = p(93);
	netValence_r12 = p(94);
	kturnover = p(95);
	Fak_init_uM = p(96);
	ksf = p(97);
	YAPTAZN_init_uM = p(98);
	n2 = p(99);
	n1 = p(100);
	NPCA_init_molecules_um_2 = p(101);
	Positionboolean_init_molecules_um_2 = p(102);
	KMOLE = p(103);
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);

    MRTFReleaseConst = p(122);
    cytoDConc = p(123);

    MRTFTot = p(124);
    Kinsolo_MRTF = p(125);
    Kin2_MRTF  = p(126);
    % kout_MRTF = p(127);
    Kdim = p(127);
    Kcap = p(128);
	
    % SS calcs
    Faktot = Fak_init_uM + Fakp_init_uM;
    Fakp = Faktot * (ksf*Emol + kf *(C+Emol)) / ((kf+kdf)*(C + Emol) + ksf*Emol);
    Fak = Faktot - Fakp;
    FakTau = (kf + kdf + ksf*Emol/(Emol + C))^(-1);

    RhoAGTP_init_uM = RhoAGTP_MEM_init_molecules_um_2 * (UnitFactor_uM_um3_molecules_neg_1*KFlux_PM_Cyto);
    RhoATot = RhoAGDP_init_uM + RhoAGTP_init_uM;
    RhoAGTP = RhoATot * kfkp*(gamma*Fakp^n2 + 1) / (kdp + kfkp*(gamma*Fakp^n2 + 1));
    RhoAGDP = RhoATot - RhoAGTP;
    RhoAGTP_MEM = RhoAGTP * (unitconversionfactor * SAV);
    RhoATau = (kdp + kfkp*(gamma*Fakp^n2 + 1))^(-1);


    ROCKTot = ROCKA_init_uM + ROCK_init_uM;
    ROCKA = ROCKTot * (krp*RhoAGTP/(kdrock + krp*RhoAGTP));
    ROCK = ROCKTot - ROCKA;
    ROCKTau = (kdrock + krp*RhoAGTP)^(-1);

    mDiaTot = mDia_init_uM + mDiaA_init_uM;
    mDiaA = mDiaTot*kmp*RhoAGTP/(kmp*RhoAGTP + kdmdia);
    mDia = mDiaTot - mDiaA;
    mDiaTau = (kmp*RhoAGTP + kdmdia)^(-1);

    % define smooth/tanh functions
    smoothROCKA = (ROCKA/2) * (tanh(20*(ROCKA-ROCKB)) + 1);
    smoothmDiaA = (mDiaA/2) * (tanh(20*(mDiaA-mDiaB)) + 1);

    LIMKTot = LIMKA_init_uM + LIMK_init_uM;
    LIMKA = LIMKTot*klr*(tau*smoothROCKA + 1)/(kdl + klr*(tau*smoothROCKA + 1));
    LIMK = LIMKTot - LIMKA;
    LIMKTau = (kdl + klr*(tau*smoothROCKA + 1))^(-1);

    CofilinTot = CofilinP_init_uM + CofilinNP_init_uM;
    CofilinRoots = roots([kturnover, ...
                       kcatcof*LIMKA + kmCof*kturnover - kturnover*CofilinTot, ...
                       -kturnover*kmCof*CofilinTot]);
    CofilinNP = CofilinRoots(CofilinRoots > 0);
    CofilinP = CofilinTot - CofilinNP;
    CofilinTau = (kturnover + kcatcof*LIMKA/kmCof)^(-1);

    ActinTot = Fcyto_init_uM + Gactin_init_uM;
    kdep_tot = (kdep + kfc1*CofilinNP)*(1 + cytoDConc/Kdim);
    Fcyto = ActinTot*kra*(alpha*smoothmDiaA+1)...
        /(kdep_tot + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/Kcap));
    Gactin = (ActinTot - Fcyto*(1+cytoDConc/Kcap))/(1 + cytoDConc/Kdim);
    ActinTau = (kdep_tot + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/Kcap))^(-1);

    MyoTot = Myo_init_uM + MyoA_init_uM;
    MyoA = MyoTot*kmr*(epsilon*smoothROCKA+1)/(kdmy + kmr*(epsilon*smoothROCKA+1));
    Myo = MyoTot - MyoA;
    MyoTau = (kdmy + kmr*(epsilon*smoothROCKA+1))^(-1);

    Ecytosol = propStiff*(Fcyto^n1);
    LaminATot = LaminAp_init_molecules_um_2 + LaminA_init_molecules_um_2;
    LaminA = LaminATot*kflaminA*Ecytosol/(kflaminA*Ecytosol + Kr_r15*(Clamin+Ecytosol));
    LaminAp = LaminATot - LaminA;
    LaminATau = (kflaminA*Ecytosol/(Clamin+Ecytosol) + Kr_r15)^(-1);

    NPCTot = NPCA_init_molecules_um_2 + NPC_init_molecules_um_2;
    NPCA = NPCTot*kfNPC*LaminA*Fcyto*MyoA / (krNPC + kfNPC*LaminA*Fcyto*MyoA);
    NPC = NPCTot - NPCA;
    NPCTau = (krNPC + kfNPC*LaminA*Fcyto*MyoA)^(-1);

    % YAPTAZ (calc in terms of molecules, then convert)
    CytoConvert = Size_Cyto/UnitFactor_uM_um3_molecules_neg_1;
    NucConvert = Size_Nuc/UnitFactor_uM_um3_molecules_neg_1;
    YAPTAZTot =  (YAPTAZP_init_uM*CytoConvert + YAPTAZN_init_uM*CytoConvert ...
        + YAPTAZnuc_init_uM*NucConvert);
    YAPTAZPFraction = (kNC/(kNC + kCN + kCY*Fcyto*MyoA)); % rapid phosphorylation
    YAPTAZnucFraction = (kinSolo2 + kin2*NPCA)/(kout2*Size_Cyto/Size_Nuc + kinSolo2 + kin2*NPCA);
%     YAPTAZCytoTot = YAPTAZTot - YAPTAZnuc;
%     YAPTAZP = YAPTAZCytoTot * YAPTAZPFraction;
%     YAPTAZN = YAPTAZCytoTot - YAPTAZP;
%     YAPTAZNPTot = YAPTAZTot - YAPTAZP;
    YAPTAZCytoTot = YAPTAZTot / (1 + (1-YAPTAZPFraction)*YAPTAZnucFraction/(1-YAPTAZnucFraction));
    YAPTAZP = YAPTAZCytoTot * YAPTAZPFraction;
    YAPTAZN = YAPTAZCytoTot - YAPTAZP;
    YAPTAZnuc = YAPTAZTot - YAPTAZCytoTot;
    if p(136) > 0 % overexpression of YAP mutant
        YAPTAZTot_5SA = YAPTAZTot*p(136);
        kNC_5SA = 0.1*kNC;
        YAPTAZPFraction_5SA = (kNC_5SA/(kNC_5SA + kCN + kCY*Fcyto*MyoA)); % rapid phosphorylation
        YAPTAZnucFraction_5SA = (kinSolo2 + kin2*NPCA)/(kout2*Size_Cyto/Size_Nuc + kinSolo2 + kin2*NPCA);
        YAPTAZCytoTot_5SA = YAPTAZTot_5SA / (1 + (1-YAPTAZPFraction_5SA)*YAPTAZnucFraction_5SA/(1-YAPTAZnucFraction_5SA));
        YAPTAZP_5SA = YAPTAZCytoTot_5SA * YAPTAZPFraction_5SA;
        YAPTAZN_5SA = YAPTAZCytoTot_5SA - YAPTAZP_5SA;
        YAPTAZnuc_5SA = YAPTAZTot_5SA - YAPTAZCytoTot_5SA;
        YAPTAZP = YAPTAZP + YAPTAZP_5SA;
        YAPTAZN = YAPTAZN + YAPTAZN_5SA;
        YAPTAZnuc = YAPTAZnuc + YAPTAZnuc_5SA;
    end

    YAPTAZN = YAPTAZN/CytoConvert;
    YAPTAZP = YAPTAZP/CytoConvert;
    YAPTAZnuc = YAPTAZnuc/NucConvert;
    YAPTAZPTau = (kNC + kCN + kCY*Fcyto*MyoA)^(-1);
    YAPTAZnucTau = (Size_Nuc*(kout2/NucConvert + (kinSolo2 + kin2*NPCA)/CytoConvert))^(-1);

    kinMod = (kinSolo2 + kin2*NPCA)*(1 - YAPTAZPFraction)*NucConvert/CytoConvert;
    YAPTAZnuc_test = (YAPTAZTot/NucConvert)*(kinMod/(kinMod+kout2));

    MRTFFactor = 1/(1+(Gactin/MRTFReleaseConst)^2);
    nucPerm_MRTF = Kinsolo_MRTF + Kin2_MRTF*NPCA;
    MRTFnuc = (MRTFTot*nucPerm_MRTF*MRTFFactor/CytoConvert) / (nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + 1);
    MRTFcyto = (MRTFTot - MRTFnuc*NucConvert)/CytoConvert;
    MRTFnucTau = 0;%(nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + kout_MRTF)^(-1);

	SSVar = [CofilinP; Fak; mDia; LaminA; Fcyto; RhoAGTP_MEM; mDiaA; NPCA;...
        Gactin; NPC; ROCKA; Myo; CofilinNP; LaminAp; YAPTAZnuc; Fakp;...
        YAPTAZP; YAPTAZN; RhoAGDP; LIMK; MyoA; ROCK; 1; LIMKA; MRTFnuc; MRTFcyto]';
    tauVals = [CofilinTau; FakTau; mDiaTau; LaminATau; ActinTau; RhoATau; mDiaTau; NPCTau;...
        ActinTau; NPCTau; ROCKTau; MyoTau; CofilinTau; LaminATau; YAPTAZnucTau; FakTau;...
        YAPTAZPTau; YAPTAZPTau; RhoATau; LIMKTau; MyoTau; ROCKTau; 1; LIMKTau; MRTFnucTau; MRTFnucTau]';
    
end