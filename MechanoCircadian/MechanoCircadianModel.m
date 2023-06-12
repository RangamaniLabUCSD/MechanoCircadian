function [T,Y, SSVals] = YAPTAZClock_EAF(timeSpan, treatmentParam, paramList, inhibVec, popVar)
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

%
% Default time span
%

% output variable lengh and names
numVars = 82;
allNames = {'CofilinP';'Fak';'mDia';'LaminA';'Fcyto';'RhoAGTP_MEM';'mDiaA';'NPCA';'Gactin';'NPC';'ROCKA';'Myo';'CofilinNP';'LaminAp';'YAPTAZnuc';'Fakp';'YAPTAZP';'YAPTAZN';'RhoAGDP';'LIMK';'MyoA';'ROCK';'Positionboolean';'LIMKA';'Kf_r9';'Kf_r8';'Kf_r7';'Kf_r6';'Kf_r5';'Kf_r4';'Kf_r3';'Kf_r2';'Kf_r1';'Kf_r0';'Ecytosol';'Kr_r9';'J_r9';'Kr_r8';'J_r8';'J_r7';'J_r6';'J_r5';'J_r4';'J_r3';'J_r2';'J_r1';'J_r0';'KFlux_NM_Nuc';'Kr_r16';'Kf_r16';'J_r16';'Kf_r15';'J_r15';'Kr_r14';'Kf_r14';'J_r14';'Kr_r13';'Kf_r13';'J_r13';'Kf_r12';'J_r12';'Kr_r11';'Kf_r11';'J_r11';'Kr_r10';'Kf_r10';'J_r10';'KFlux_PM_Cyto';'KFlux_NM_Cyto';'UnitFactor_uM_um3_molecules_neg_1';'YAPTAZratio';'FakConservation';'RhoAConservation';'LIMKConservation';'YAPTAZConservation';'MyoConservation';'mDiaConservation';'ROCKConservation';'CofilinConservation';'LaminAConservation';'NPCConservation';'ActinConservation'};

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
param(17) = treatmentParam(1); % substrate stiffness
param(104) = treatmentParam(2); % time parameter for stiffness increase

varVecCur = popVar.*randn(size(param(1:104)));
varVecCur([33:36,47:49,51,64,66,68,72,85,90:92,94,99:100,102:103]) = 0;
param(1:104) = param(1:104) .* exp(varVecCur);

if length(paramList) ~= 19 % then default parameters
    param(105) = 9.56*3600; % tauB1
    param(106) = 3.29; % nExpB
    param(107) = 1.29/3600; % KeB
    param(108) = 0.05; % KiB
    param(109) = 0.62/3600; % KdBP
    param(110) = 0.49/3600; % KdB
    param(111) = 11.05*3600; % tauB2
    param(112) = 3.04; % nExpP
    param(113) = 1.29/3600; % KeP
    param(114) = 0.63; % KaP
    param(115) = 0.50/3600; % KdP

    param(116) = 0.05/3600; %magn of YAPTAZ->Ke coupling
    param(117) = 1; %Kd for YAPTAZ->Ke
    param(118) = 2; % Hill coeff for YAPTAZ->Ke
    param(119) = 0.05/3600;
    param(120) = 1;
    param(121) = 2;
    param(122) = 450;
else
    param(105:121) = paramList(1:17);
    % load inhibition parameters [actin polym (kra), ROCK, MRTF, YAP phosphorylation (kNC)]
    param(3) = param(3)*inhibVec(1); % inhibition of actin polymerization
    param([55,69]) = param([55,69])*inhibVec(2); % inhibit ROCK-mediated catalysis
    param(119) = param(119)*inhibVec(3); % inhibit coupling between MRTF and circadian
    param(2) = param(2)*inhibVec(4); % inhibit YAP phosphorylation
    if inhibVec(4) < 1
        param([4,56,98]) = 2*param([4,56,98]);
    end
    % if applicable, assess cytoD conc
    param(123) = 0;
    if length(inhibVec)>4
        param(123) = inhibVec(5);
    end
    if length(inhibVec)>5
        param(2) = inhibVec(6)*param(2);
    end
    param(78) = paramList(19);
end

param(122) = 43.2495; %112.5; %MRTFReleaseConst
param(124) = 1.0724e7; %1.0153; %MRTFTot
param(125) = 8.4028; %2.6383; %kinsolo_MRTF
param(126) = 16.5570; %6.2886; %kin2_MRTF
param(127) = 0.1346; %0.2725; %kout_MRTF
param(128) = 0.3383; %0.25; %cytoDConst


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
[T,Y] = fDDE(timeSpan,param);
SSVals = fSS(param);


% % get the solution
% all = zeros(length(T), numVars);
% for i = 1:length(T)
% 	all(i,:) = getRow(T(i), Y(i,:), yinit, param);
% end
% 
% allValues = all;
end

% -------------------------------------------------------
% get row data
function rowValue = getRow(t,y,y0,p)
	% State Variables
	CofilinP = y(1);
	Fak = y(2);
	mDia = y(3);
	LaminA = y(4);
	Fcyto = y(5);
	RhoAGTP_MEM = y(6);
	mDiaA = y(7);
	NPCA = y(8);
	Gactin = y(9);
	NPC = y(10);
	ROCKA = y(11);
	Myo = y(12);
	CofilinNP = y(13);
	LaminAp = y(14);
	YAPTAZnuc = y(15);
	Fakp = y(16);
	YAPTAZP = y(17);
	YAPTAZN = y(18);
	RhoAGDP = y(19);
	LIMK = y(20);
	MyoA = y(21);
	ROCK = y(22);
	Positionboolean = y(23);
	LIMKA = y(24);
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
	% Functions
	Kf_r9 = (kmr .* (1.0 + (0.5 .* epsilon .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	Kf_r8 = (klr .* (1.0 + (0.5 .* tau .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	Kf_r7 = kdmdia;
	Kf_r6 = kdrock;
	Kf_r5 = (kmp .* RhoAGTP_MEM);
	Kf_r4 = (krp .* RhoAGTP_MEM);
	Kf_r3 = kdp;
	Kf_r2 = kf;
	Kf_r1 = (kfkp .* (1.0 + (gamma .* (Fakp ^ n2))) .* unitconversionfactor .* SAV);
	Kf_r0 = kdf;
	Ecytosol = (propStiff .* (Fcyto ^ n1));
	Kr_r9 = kdmy;
	J_r9 = ((Kf_r9 .* Myo) - (Kr_r9 .* MyoA));
	Kr_r8 = kdl;
	J_r8 = ((Kf_r8 .* LIMK) - (Kr_r8 .* LIMKA));
	J_r7 = ((Kf_r7 .* mDiaA) - (Kr_r7 .* mDia));
	J_r6 = ((Kf_r6 .* ROCKA) - (Kr_r6 .* ROCK));
	J_r5 = ((Kf_r5 .* mDia) - (Kr_r5 .* mDiaA));
	J_r4 = ((Kf_r4 .* ROCK) - (Kr_r4 .* ROCKA));
	J_r3 = ((Kf_r3 .* RhoAGTP_MEM) - (Kr_r3 .* RhoAGDP));
	J_r2 = ((Kf_r2 .* Fak) - (Kr_r2 .* Fakp));
	J_r1 = ((Kf_r1 .* RhoAGDP) - (Kr_r1 .* RhoAGTP_MEM));
	J_r0 = ((Kf_r0 .* Fakp) - (Kr_r0 .* Fak));
	KFlux_NM_Nuc = (Size_NM ./ Size_Nuc);
	Kr_r16 = kout2;
	Kf_r16 = ((kin2 .* NPCA) + kinSolo2);
	J_r16 = ((Kf_r16 .* YAPTAZN) - (Kr_r16 .* YAPTAZnuc));
	Kf_r15 = (kflaminA .* Ecytosol ./ (Clamin + Ecytosol));
	J_r15 = ((Kf_r15 .* LaminAp) - (Kr_r15 .* LaminA));
	Kr_r14 = krNPC;
	Kf_r14 = (kfNPC .* MyoA .* Fcyto .* LaminA);
	J_r14 = ((Kf_r14 .* NPC) - (Kr_r14 .* NPCA));
	Kr_r13 = kNC;
	Kf_r13 = (kCN + (kCY .* MyoA .* Fcyto));
	J_r13 = ((Kf_r13 .* YAPTAZP) - (Kr_r13 .* YAPTAZN));
	Kf_r12 = (ksf .* Emol .* unitconversionfactor .* SAV .* Positionboolean ./ (C + Emol));
	J_r12 = ((Kf_r12 .* Fak) - (Kr_r12 .* Fakp));
	Kr_r11 = (kdep + (kfc1 .* CofilinNP));
	Kf_r11 = (kra .* (1.0 + (0.5 .* alpha .* (1.0 + tanh((20.0 .* (mDiaA - mDiaB)))) .* mDiaA)));
	J_r11 = ((Kf_r11 .* Gactin) - (Kr_r11 .* Fcyto));
	Kr_r10 = (kcatcof .* LIMKA ./ (kmCof + CofilinNP));
	Kf_r10 = kturnover;
	J_r10 = ((Kf_r10 .* CofilinP) - (Kr_r10 .* CofilinNP));
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);
	% OutputFunctions
	YAPTAZratio = (YAPTAZnuc ./ (YAPTAZN + YAPTAZP));
	FakConservation = (Fak + Fakp);
	RhoAConservation = (RhoAGDP + (RhoAGTP_MEM .* 1256.0 ./ 2300.0 .* 1.0E15 ./ 6.022E17));
	LIMKConservation = (LIMKA + LIMK);
	YAPTAZConservation = (((YAPTAZN + YAPTAZP) .* 2300.0) + (YAPTAZnuc .* 550.0));
	MyoConservation = (Myo + MyoA);
	mDiaConservation = (mDia + mDiaA);
	ROCKConservation = (ROCK + ROCKA);
	CofilinConservation = (CofilinNP + CofilinP);
	LaminAConservation = (LaminAp + LaminA);
	NPCConservation = (NPC + NPCA);
	ActinConservation = (Fcyto + Gactin);

	rowValue = [CofilinP Fak mDia LaminA Fcyto RhoAGTP_MEM mDiaA NPCA Gactin NPC ROCKA Myo CofilinNP LaminAp YAPTAZnuc Fakp YAPTAZP YAPTAZN RhoAGDP LIMK MyoA ROCK Positionboolean LIMKA Kf_r9 Kf_r8 Kf_r7 Kf_r6 Kf_r5 Kf_r4 Kf_r3 Kf_r2 Kf_r1 Kf_r0 Ecytosol Kr_r9 J_r9 Kr_r8 J_r8 J_r7 J_r6 J_r5 J_r4 J_r3 J_r2 J_r1 J_r0 KFlux_NM_Nuc Kr_r16 Kf_r16 J_r16 Kf_r15 J_r15 Kr_r14 Kf_r14 J_r14 Kr_r13 Kf_r13 J_r13 Kf_r12 J_r12 Kr_r11 Kf_r11 J_r11 Kr_r10 Kf_r10 J_r10 KFlux_PM_Cyto KFlux_NM_Cyto UnitFactor_uM_um3_molecules_neg_1 YAPTAZratio FakConservation RhoAConservation LIMKConservation YAPTAZConservation MyoConservation mDiaConservation ROCKConservation CofilinConservation LaminAConservation NPCConservation ActinConservation];
end

% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0)
	% State Variables
	CofilinP = y(1);
	Fak = y(2);
	mDia = y(3);
	LaminA = y(4);
	Fcyto = y(5);
	RhoAGTP_MEM = y(6);
	mDiaA = y(7);
	NPCA = y(8);
	Gactin = y(9);
	NPC = y(10);
	ROCKA = y(11);
	Myo = y(12);
	CofilinNP = y(13);
	LaminAp = y(14);
	YAPTAZnuc = y(15);
	Fakp = y(16);
	YAPTAZP = y(17);
	YAPTAZN = y(18);
	RhoAGDP = y(19);
	LIMK = y(20);
	MyoA = y(21);
	ROCK = y(22);
	Positionboolean = y(23);
	LIMKA = y(24);
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
    if isinf(p(104))
        Emol = p(17);
    else
        Emol = p(17)*t/p(104);
    end
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
	% Functions
	Kf_r9 = (kmr .* (1.0 + (0.5 .* epsilon .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	Kf_r8 = (klr .* (1.0 + (0.5 .* tau .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	Kf_r7 = kdmdia;
	Kf_r6 = kdrock;
	Kf_r5 = (kmp .* RhoAGTP_MEM);
	Kf_r4 = (krp .* RhoAGTP_MEM);
	Kf_r3 = kdp;
	Kf_r2 = kf;
	Kf_r1 = (kfkp .* (1.0 + (gamma .* (Fakp ^ n2))) .* unitconversionfactor .* SAV);
	Kf_r0 = kdf;
	Ecytosol = (propStiff .* (Fcyto ^ n1));
	Kr_r9 = kdmy;
	J_r9 = ((Kf_r9 .* Myo) - (Kr_r9 .* MyoA));
	Kr_r8 = kdl;
	J_r8 = ((Kf_r8 .* LIMK) - (Kr_r8 .* LIMKA));
	J_r7 = ((Kf_r7 .* mDiaA) - (Kr_r7 .* mDia));
	J_r6 = ((Kf_r6 .* ROCKA) - (Kr_r6 .* ROCK));
	J_r5 = ((Kf_r5 .* mDia) - (Kr_r5 .* mDiaA));
	J_r4 = ((Kf_r4 .* ROCK) - (Kr_r4 .* ROCKA));
	J_r3 = ((Kf_r3 .* RhoAGTP_MEM) - (Kr_r3 .* RhoAGDP));
	J_r2 = ((Kf_r2 .* Fak) - (Kr_r2 .* Fakp));
	J_r1 = ((Kf_r1 .* RhoAGDP) - (Kr_r1 .* RhoAGTP_MEM));
	J_r0 = ((Kf_r0 .* Fakp) - (Kr_r0 .* Fak));
	KFlux_NM_Nuc = (Size_NM ./ Size_Nuc);
	Kr_r16 = kout2;
	Kf_r16 = ((kin2 .* NPCA) + kinSolo2);
	J_r16 = ((Kf_r16 .* YAPTAZN) - (Kr_r16 .* YAPTAZnuc));
	Kf_r15 = (kflaminA .* Ecytosol ./ (Clamin + Ecytosol));
	J_r15 = ((Kf_r15 .* LaminAp) - (Kr_r15 .* LaminA));
	Kr_r14 = krNPC;
	Kf_r14 = (kfNPC .* MyoA .* Fcyto .* LaminA);
	J_r14 = ((Kf_r14 .* NPC) - (Kr_r14 .* NPCA));
	Kr_r13 = kNC;
	Kf_r13 = (kCN + (kCY .* MyoA .* Fcyto));
	J_r13 = ((Kf_r13 .* YAPTAZP) - (Kr_r13 .* YAPTAZN));
	Kf_r12 = (ksf .* Emol .* unitconversionfactor .* SAV .* Positionboolean ./ (C + Emol));
	J_r12 = ((Kf_r12 .* Fak) - (Kr_r12 .* Fakp));
	Kr_r11 = (kdep + (kfc1 .* CofilinNP));
	Kf_r11 = (kra .* (1.0 + (0.5 .* alpha .* (1.0 + tanh((20.0 .* (mDiaA - mDiaB)))) .* mDiaA)));
	J_r11 = ((Kf_r11 .* Gactin) - (Kr_r11 .* Fcyto));
	Kr_r10 = (kcatcof .* LIMKA ./ (kmCof + CofilinNP));
	Kf_r10 = kturnover;
	J_r10 = ((Kf_r10 .* CofilinP) - (Kr_r10 .* CofilinNP));
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);
	% Rates
	dydt = [
		 - J_r10;    % rate for CofilinP
		( - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r12) + J_r0 - J_r2);    % rate for Fak
		( - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r5) + J_r7);    % rate for mDia
		J_r15;    % rate for LaminA
		J_r11;    % rate for Fcyto
		( - J_r3 + J_r1);    % rate for RhoAGTP_MEM
		((UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r5) - J_r7);    % rate for mDiaA
		J_r14;    % rate for NPCA
		 - J_r11;    % rate for Gactin
		 - J_r14;    % rate for NPC
		((UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r4) - J_r6);    % rate for ROCKA
		 - J_r9;    % rate for Myo
		J_r10;    % rate for CofilinNP
		 - J_r15;    % rate for LaminAp
		(UnitFactor_uM_um3_molecules_neg_1 .* KFlux_NM_Nuc .* J_r16);    % rate for YAPTAZnuc
		((UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r12) - J_r0 + J_r2);    % rate for Fakp
		 - J_r13;    % rate for YAPTAZP
		(J_r13 - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_NM_Cyto .* J_r16));    % rate for YAPTAZN
		((UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r3) - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r1));    % rate for RhoAGDP
		 - J_r8;    % rate for LIMK
		J_r9;    % rate for MyoA
		( - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_PM_Cyto .* J_r4) + J_r6);    % rate for ROCK
		0.0;    % rate for Positionboolean
		J_r8;    % rate for LIMKA
	];
end


% reduced dde
function [t,y] = fDDE(timeSpan,p)
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
    tauB1 = p(105);
    pExpB = p(106);
    KeB = p(107);
    KiB = p(108);
    KdBP = p(109);
    KdB = p(110);
    tauB2 = p(111);
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
    kinsolo_MRTF = p(125);
    kin2_MRTF  = p(126);
    kout_MRTF = p(127);
    cytoDConst = p(128);
    
    lags = [tauB1, tauB1];
%     lags = [tauB1, tauB2];

%     DDESol = dde23(@ddefun, lags, @history, timeSpan);
    SSVals = fSS(p);
    YAPTAZnuc_SS = SSVals(15);
    MRTFnuc_SS = SSVals(25);
    DDESol = dde23(@ddefun_YAPSS, lags, @history_CircOnly, timeSpan);
    t = DDESol.x';
    y = DDESol.y';
	
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
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc)^nCouple1);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc)^nCouple2);

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
        KeB2 = magCouple1 / (1 + (KdCouple1/YAPTAZnuc_SS)^nCouple1);
        KeP2 = magCouple2 / (1 + (KdCouple2/MRTFnuc_SS)^nCouple2);

        % Rates (BMAL and PER dynamics)
        dydt = [KeB*(1./(1+(clockLag1/KiB).^pExpB)) - KdBP*y(1)*y(2) - KdB*y(1) + KeB2;
            KeP*(1./(1+(KaP/clockLag2).^pExpP)) - KdBP*y(1)*y(2) - KdP*y(2) + KeP2];
    end

    function s = history_CircOnly(t)
        s = [0;
             0];
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
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);

    MRTFReleaseConst = p(122);
    cytoDConc = p(123);

    MRTFTot = p(124);
    kinsolo_MRTF = p(125);
    kin2_MRTF  = p(126);
    kout_MRTF = p(127);
    cytoDConst = p(128);
	
    %SS calcs (some analytical)
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
    Fcyto = ActinTot*kra*(alpha*smoothmDiaA+1)...
        /(kdep + kfc1*CofilinNP + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/cytoDConst));
    Gactin = ActinTot - Fcyto*(1+cytoDConc/cytoDConst);
    ActinTau = (kdep + kfc1*CofilinNP + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/cytoDConst))^(-1);

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
    YAPTAZN = YAPTAZN/CytoConvert;
    YAPTAZP = YAPTAZP/CytoConvert;
    YAPTAZnuc = YAPTAZnuc/NucConvert;
    YAPTAZPTau = (kNC + kCN + kCY*Fcyto*MyoA)^(-1);
    YAPTAZnucTau = (Size_Nuc*(kout2/NucConvert + (kinSolo2 + kin2*NPCA)/CytoConvert))^(-1);

    MRTFFactor = 1/(1+(Gactin/MRTFReleaseConst)^2);
    nucPerm_MRTF = kinsolo_MRTF + kin2_MRTF*NPCA;
    MRTFnuc = (MRTFTot*nucPerm_MRTF*MRTFFactor/CytoConvert) / (nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + kout_MRTF);
    MRTFcyto = (MRTFTot - MRTFnuc*NucConvert)/CytoConvert;
    MRTFnucTau = (nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + kout_MRTF)^(-1);

	SSVar = [CofilinP; Fak; mDia; LaminA; Fcyto; RhoAGTP_MEM; mDiaA; NPCA;...
        Gactin; NPC; ROCKA; Myo; CofilinNP; LaminAp; YAPTAZnuc; Fakp;...
        YAPTAZP; YAPTAZN; RhoAGDP; LIMK; MyoA; ROCK; 1; LIMKA; MRTFnuc; MRTFcyto]';
    tauVals = [CofilinTau; FakTau; mDiaTau; LaminATau; ActinTau; RhoATau; mDiaTau; NPCTau;...
        ActinTau; NPCTau; ROCKTau; MyoTau; CofilinTau; LaminATau; YAPTAZnucTau; FakTau;...
        YAPTAZPTau; YAPTAZPTau; RhoATau; LIMKTau; MyoTau; ROCKTau; 1; LIMKTau; MRTFnucTau; MRTFnucTau]';
    
end