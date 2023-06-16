function [T,Y,yinit,param, allNames, allValues] = MechanoCircadianModel_full(argTimeSpan,argYinit,argParam,stiffnessVal,couplingParam)
% [T,Y,yinit,param] = YAPTAZ_FullClockNew(argTimeSpan,argYinit,argParam)
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
timeSpan = [0.0 1.0];

% output variable lengh and names
numVars = 254;
allNames = {'MRTF';'PERCRYnuc';'LaminA';'MyoA';'BMAL1';'Fakp';'TEAD';'PERCRY_CLOCKBMAL1';'Gactin';'Fak';'ROCK';'BMAL1_mRNA';'CofilinNP';'REV';'BMAL1nuc';'NPCA';'PER_mRNA';'Positionboolean';'LIMK';'RhoAGDP';'MRTFnuc';'Fcyto';'BMAL1nuc_p';'Myo';'PERCRY_p';'YAPTAZP';'YAPTAZN';'REV_mRNA';'SRF';'CRY';'LaminAp';'PERCRYnuc_p';'YAPTAZnuc';'PER_p';'mDiaA';'LIMKA';'mDia';'PERCRY';'CRY_mRNA';'ROCKA';'RhoAGTP_MEM';'PER';'REVnuc';'CofilinP';'CRY_p';'NPC';'BMAL1_p';'UnitFactor_uM_um3_molecules_neg_1';'nu_dRC';'K_d_r36';'K_d_r32';'k_sC';'LumpedJ_r44';'K_d_r30';'k_sP';'LumpedJ_r42';'KFlux_PM_Cyto';'k_9';'k_6';'k_5';'k_dnC';'J_r39';'K_d_r38';'k_dn_r38';'nu_dCC';'J_r38';'k_dn_r37';'J_r37';'k_2';'nu_dPCC';'k_dn_r36';'J_r36';'k_1';'V_1P';'K_dp_r35';'K_p_r35';'V_2P';'J_r35';'V_2C';'V_1C';'K_p_r34';'K_dp_r34';'J_r34';'K_dp_r33';'K_p_r33';'V_2PC';'V_1PC';'J_r33';'nu_dBN';'k_dn_r32';'J_r32';'k_dn_r31';'J_r31';'k_dn_r30';'nu_dBC';'J_r30';'K_d_r26';'K_d_r24';'Kr_r9';'Kf_r9';'J_r9';'Kr_r8';'Kf_r8';'J_r8';'Kf_r7';'J_r7';'Kf_r6';'J_r6';'Kf_r5';'J_r5';'Kf_r4';'J_r4';'Kf_r3';'J_r3';'k_dn_r29';'J_r29';'k_dmB';'nu_mB';'K_mB';'J_r28';'Kf_r2';'J_r2';'nu_sB';'nu_YTB';'K_IB';'J_r27';'Kf_r1';'J_r1';'k_dn_r26';'nu_dRN';'J_r26';'Kf_r0';'J_r0';'k_dmr';'nu_mR';'K_mR';'J_r25';'k_dn_r24';'J_r24';'nu_sR';'K_AR';'J_r23';'V_2B';'V_1B';'K_dp_r20';'K_p_r20';'J_r20';'nu_dPC';'k_10';'Kr_r22';'Kf_r22';'LumpedJ_r22';'LumpedJ_r21';'K_p_r52';'V_3B';'K_p_r19';'K_dp_r19';'V_4B';'J_r19';'KFlux_NM_Nuc';'Kr_r16';'Kf_r16';'J_r16';'Ecytosol';'Kf_r15';'J_r15';'Kr_r14';'Kf_r14';'J_r14';'Kr_r13';'Kf_r13';'J_r13';'Kf_r12';'J_r12';'Kr_r11';'Kf_r11';'J_r11';'Kr_r10';'Kf_r10';'J_r10';'k_sB';'LumpedJ_r18';'k_sR';'LumpedJ_r17';'K_AP';'K_AC';'nu_MRTF_r49';'nu_MRTF_r46';'k_4_r51';'nu_sP';'k_4_r45';'k_4_r43';'nu_sC';'k_3_r51';'nu_dPCN';'k_3_r45';'k_3_r43';'nu_dIN';'K_mP';'K_mC';'nu_mP';'nu_mC';'KFlux_NM_Cyto';'k_dn_r58';'k_dn_r54';'k_dn_r53';'k_dn_r41';'k_dn_r40';'K_dp_r52';'Kf_r59';'V_4PC';'Kr_r59';'k_dmp';'k_dmc';'k_8_r57';'K_d_r58';'k_8_r56';'k_8_r55';'K_d_r54';'J_r59';'V_3PC';'J_r58';'k_7_r57';'J_r57';'k_7_r56';'J_r56';'k_7_r55';'J_r55';'J_r54';'J_r53';'J_r52';'J_r51';'K_d_r40';'LumpedJ_r50';'J_r49';'J_r48';'J_r47';'J_r46';'J_r45';'J_r43';'J_r41';'J_r40';'YAPTAZratio';'FakConservation';'RhoAConservation';'LIMKConservation';'YAPTAZConservation';'MyoConservation';'mDiaConservation';'ROCKConservation';'CofilinConservation';'LaminAConservation';'NPCConservation';'ActinConservation'};

if nargin >= 1
	if length(argTimeSpan) > 0
		%
		% TimeSpan overridden by function arguments
		%
		timeSpan = argTimeSpan;
	end
end
%
% Default Initial Conditions
%
yinit = [
	1.0;		% yinit(1) is the initial condition for 'MRTF'
	0.002;		% yinit(2) is the initial condition for 'PERCRYnuc'
	0.0;		% yinit(3) is the initial condition for 'LaminA'
	1.5;		% yinit(4) is the initial condition for 'MyoA'
	0.005;		% yinit(5) is the initial condition for 'BMAL1'
	0.3;		% yinit(6) is the initial condition for 'Fakp'
	0.0;		% yinit(7) is the initial condition for 'TEAD'
	0.0;		% yinit(8) is the initial condition for 'PERCRY_CLOCKBMAL1'
	482.4;		% yinit(9) is the initial condition for 'Gactin'
	0.7;		% yinit(10) is the initial condition for 'Fak'
	1.0;		% yinit(11) is the initial condition for 'ROCK'
	0.002;		% yinit(12) is the initial condition for 'BMAL1_mRNA'
	1.8;		% yinit(13) is the initial condition for 'CofilinNP'
	0.0;		% yinit(14) is the initial condition for 'REV'
	0.005;		% yinit(15) is the initial condition for 'BMAL1nuc'
	0.0;		% yinit(16) is the initial condition for 'NPCA'
	0.002;		% yinit(17) is the initial condition for 'PER_mRNA'
	1.0;		% yinit(18) is the initial condition for 'Positionboolean'
	1.9;		% yinit(19) is the initial condition for 'LIMK'
	1.0;		% yinit(20) is the initial condition for 'RhoAGDP'
	0.0;		% yinit(21) is the initial condition for 'MRTFnuc'
	17.9;		% yinit(22) is the initial condition for 'Fcyto'
	0.0;		% yinit(23) is the initial condition for 'BMAL1nuc_p'
	3.5;		% yinit(24) is the initial condition for 'Myo'
	0.0;		% yinit(25) is the initial condition for 'PERCRY_p'
	0.2;		% yinit(26) is the initial condition for 'YAPTAZP'
	0.7;		% yinit(27) is the initial condition for 'YAPTAZN'
	0.002;		% yinit(28) is the initial condition for 'REV_mRNA'
	0.0;		% yinit(29) is the initial condition for 'SRF'
	0.002;		% yinit(30) is the initial condition for 'CRY'
	3500.0;		% yinit(31) is the initial condition for 'LaminAp'
	0.0;		% yinit(32) is the initial condition for 'PERCRYnuc_p'
	0.7;		% yinit(33) is the initial condition for 'YAPTAZnuc'
	0.0;		% yinit(34) is the initial condition for 'PER_p'
	0.0;		% yinit(35) is the initial condition for 'mDiaA'
	0.1;		% yinit(36) is the initial condition for 'LIMKA'
	0.8;		% yinit(37) is the initial condition for 'mDia'
	0.002;		% yinit(38) is the initial condition for 'PERCRY'
	0.002;		% yinit(39) is the initial condition for 'CRY_mRNA'
	0.0;		% yinit(40) is the initial condition for 'ROCKA'
	33.6;		% yinit(41) is the initial condition for 'RhoAGTP_MEM'
	0.002;		% yinit(42) is the initial condition for 'PER'
	0.0;		% yinit(43) is the initial condition for 'REVnuc'
	0.2;		% yinit(44) is the initial condition for 'CofilinP'
	0.0;		% yinit(45) is the initial condition for 'CRY_p'
	6.5;		% yinit(46) is the initial condition for 'NPC'
	0.0;		% yinit(47) is the initial condition for 'BMAL1_p'
];
if nargin >= 2
	if length(argYinit) > 0
		%
		% initial conditions overridden by function arguments
		%
		yinit = argYinit;
	end
end
%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%
param = [
	1.0;		% param(1) is 'K_MRTF_r46'
	0.0;		% param(2) is 'Voltage_NM'
	96485.3321;		% param(3) is 'mlabfix_F_'
	0.002;		% param(4) is 'CRY_init_uM'
	0.7;		% param(5) is 'YAPTAZN_init_uM'
	0.04;		% param(6) is 'kturnover'
	0.002;		% param(7) is 'PER_init_uM'
	17.9;		% param(8) is 'Fcyto_init_uM'
	1.0;		% param(9) is 'netValence_r59'
	1.0;		% param(10) is 'netValence_r50'
	0.0;		% param(11) is 'CRY_p_init_uM'
	0.8;		% param(12) is 'kdrock'
	0.002;		% param(13) is 'PERCRY_init_uM'
	0.0;		% param(14) is 'BMAL1_p_init_uM'
	0.0;		% param(15) is 'ROCKA_init_uM'
	7.6E-4;		% param(16) is 'kCY'
	1.0;		% param(17) is 'netValence_r44'
	1.0;		% param(18) is 'netValence_r42'
	0.56;		% param(19) is 'kCN'
	1.0;		% param(20) is 'Positionboolean_init_molecules_um_2'
	482.4;		% param(21) is 'Gactin_init_uM'
	0.379;		% param(22) is 'ksf'
	4.0;		% param(23) is 'kfc1'
	0.067;		% param(24) is 'kdmy'
	0.648;		% param(25) is 'krp'
	0.0;		% param(26) is 'PER_p_init_uM'
	0.4;		% param(27) is 'kra'
	1.0;		% param(28) is 'netValence_r22'
	1.0;		% param(29) is 'netValence_r21'
	77.56;		% param(30) is 'gamma'
	1.8;		% param(31) is 'SAVol'
	1.0;		% param(32) is 'netValence_r18'
	1.0;		% param(33) is 'netValence_r17'
	1.0;		% param(34) is 'netValence_r16'
	1.0;		% param(35) is 'netValence_r15'
	1.0;		% param(36) is 'netValence_r14'
	0.0;		% param(37) is 'Kr_r7'
	0.0;		% param(38) is 'Kr_r6'
	1.0;		% param(39) is 'netValence_r12'
	0.0;		% param(40) is 'Kr_r5'
	0.34;		% param(41) is 'kcatcof'
	0.0;		% param(42) is 'Kr_r4'
	0.0;		% param(43) is 'Kr_r3'
	0.002;		% param(44) is 'PERCRYnuc_init_uM'
	0.0;		% param(45) is 'Kr_r2'
	0.0;		% param(46) is 'Kr_r1'
	0.0;		% param(47) is 'Kr_r0'
	1.8;		% param(48) is 'SAV'
	1000.0;		% param(49) is 'K_millivolts_per_volt'
	1.0E-9;		% param(50) is 'mlabfix_K_GHK_'
	1260.0;		% param(51) is 'Size_PM'
	1.0;		% param(52) is 'q_r49'
	2.0;		% param(53) is 'q_r46'
	1.5;		% param(54) is 'MyoA_init_uM'
	8.4891;		% param(55) is 'kinsolo_MRTF'
	1.0;		% param(56) is 'netValence_r5'
	1.0;		% param(57) is 'netValence_r4'
	1.0;		% param(58) is 'netValence_r3'
	1.0;		% param(59) is 'netValence_r1'
	0.0;		% param(60) is 'MRTFnuc_init_uM'
	9.64853321E-5;		% param(61) is 'mlabfix_F_nmol_'
	0.0;		% param(62) is 'NPCA_init_molecules_um_2'
	393.0;		% param(63) is 'Size_NM'
	0.0;		% param(64) is 'TEAD_init_uM'
	2.0;		% param(65) is 'q_r27'
	4.0;		% param(66) is 'kmCof'
	0.03;		% param(67) is 'kmr'
	0.002;		% param(68) is 'kmp'
	0.0;		% param(69) is 'REV_init_uM'
	0.8;		% param(70) is 'mDia_init_uM'
	120.0;		% param(71) is 'kly'
	43.195;		% param(72) is 'MRTFReleaseConst'
	0.07;		% param(73) is 'klr'
	1.8;		% param(74) is 'CofilinNP_init_uM'
	16.0;		% param(75) is 'kll'
	0.1356;		% param(76) is 'kout_MRTF'
	300.0;		% param(77) is 'mlabfix_T_'
	2300.0;		% param(78) is 'Size_Cyto'
	9.0E-6;		% param(79) is 'propStiff'
	550.0;		% param(80) is 'Size_Nuc'
	100.0;		% param(81) is 'Clamin'
	5.0;		% param(82) is 'n2'
	0.0;		% param(83) is 'LaminA_init_molecules_um_2'
	2.6;		% param(84) is 'n1'
	3.5;		% param(85) is 'kdep'
	2.0;		% param(86) is 'n_r49'
	2.0;		% param(87) is 'n_r46'
	33.6;		% param(88) is 'RhoAGTP_MEM_init_molecules_um_2'
	0.002;		% param(89) is 'BMAL1_mRNA_init_uM'
	8314.46261815;		% param(90) is 'mlabfix_R_'
	0.015;		% param(91) is 'kf'
	0.2;		% param(92) is 'CofilinP_init_uM'
	0.0;		% param(93) is 'PERCRY_p_init_uM'
	1.9;		% param(94) is 'LIMK_init_uM'
	0.7;		% param(95) is 'YAPTAZnuc_init_uM'
	3500.0;		% param(96) is 'LaminAp_init_molecules_um_2'
	0.002;		% param(97) is 'REV_mRNA_init_uM'
	0.46;		% param(98) is 'kflaminA'
	55.49;		% param(99) is 'tau'
	0.0;		% param(100) is 'mDiaA_init_uM'
	0.0;		% param(101) is 'BMAL1nuc_p_init_uM'
	0.002;		% param(102) is 'CRY_mRNA_init_uM'
	0.005;		% param(103) is 'BMAL1_init_uM'
	8.7;		% param(104) is 'krNPC'
	0.0;		% param(105) is 'PERCRYnuc_p_init_uM'
	0.3;		% param(106) is 'Fakp_init_uM'
	0.14;		% param(107) is 'kNC'
	1.0;		% param(108) is 'RhoAGDP_init_uM'
	1.0;		% param(109) is 'convertBool'
	0.165;		% param(110) is 'mDiaB'
	0.0;		% param(111) is 'REVnuc_init_uM'
	0.7;		% param(112) is 'Fak_init_uM'
	3.5;		% param(113) is 'Myo_init_uM'
	0.625;		% param(114) is 'kdp'
	0.0168;		% param(115) is 'kfkp'
	2.0;		% param(116) is 'kdl'
	0.035;		% param(117) is 'kdf'
	6.5;		% param(118) is 'NPC_init_molecules_um_2'
	0.1;		% param(119) is 'LIMKA_init_uM'
	stiffnessVal;		% param(120) is 'Emol'
	2.8E-7;		% param(121) is 'kfNPC'
	3.141592653589793;		% param(122) is 'mlabfix_PI_'
	0.0;		% param(123) is 'SRF_init_uM'
	1704.7954979572205;		% param(124) is 'Size_ECM'
	50.0;		% param(125) is 'alpha'
	0.3;		% param(126) is 'ROCKB'
	0.002;		% param(127) is 'PER_mRNA_init_uM'
	602.2;		% param(128) is 'unitconversionfactor'
	2.0;		% param(129) is 'm'
	0.005;		% param(130) is 'kdmdia'
	2.0;		% param(131) is 'h'
	0.001660538783162726;		% param(132) is 'KMOLE'
	16.7249;		% param(133) is 'kin2_MRTF'
	6.02214179E11;		% param(134) is 'mlabfix_N_pmol_'
	0.2;		% param(135) is 'YAPTAZP_init_uM'
	1.0;		% param(136) is 'kinSolo2'
	16.0;		% param(137) is 'k11'
	0.1;		% param(138) is 'MRTF_init_uM'
	3.25;		% param(139) is 'C'
	1.0;		% param(140) is 'kout2'
	0.0;		% param(141) is 'Voltage_PM'
	2.0;		% param(142) is 'K_YTB'
	36.0;		% param(143) is 'epsilon'
	10.0;		% param(144) is 'kin2'
	0.005;		% param(145) is 'BMAL1nuc_init_uM'
	1.0;		% param(146) is 'ROCK_init_uM'
	0.0;		% param(147) is 'PERCRY_CLOCKBMAL1_init_uM'
	0.001;		% param(148) is 'Kr_r15'
	1.0;		% param(149) is 'K_MRTF_r49'
	0.0;		% param(150) is 'Kr_r12'
];
if nargin >= 3
	if length(argParam) > 0
		%
		% parameter values overridden by function arguments
		%
		param = argParam;
	end
end
if nargin <= 3
    stiffnessVal = 25;
end
    
%
% invoke the integrator
%
options = odeset;
[T,Y] = ode15s(@f,timeSpan,yinit,options,param,yinit,couplingParam);

% get the solution
all = zeros(length(T), numVars);
for i = 1:length(T)
	all(i,:) = getRow(T(i), Y(i,:), yinit, param);
end

allValues = all;
end

% -------------------------------------------------------
% get row data
function rowValue = getRow(t,y,y0,p)
	% State Variables
	MRTF = y(1);
	PERCRYnuc = y(2);
	LaminA = y(3);
	MyoA = y(4);
	BMAL1 = y(5);
	Fakp = y(6);
	TEAD = y(7);
	PERCRY_CLOCKBMAL1 = y(8);
	Gactin = y(9);
	Fak = y(10);
	ROCK = y(11);
	BMAL1_mRNA = y(12);
	CofilinNP = y(13);
	REV = y(14);
	BMAL1nuc = y(15);
	NPCA = y(16);
	PER_mRNA = y(17);
	Positionboolean = y(18);
	LIMK = y(19);
	RhoAGDP = y(20);
	MRTFnuc = y(21);
	Fcyto = y(22);
	BMAL1nuc_p = y(23);
	Myo = y(24);
	PERCRY_p = y(25);
	YAPTAZP = y(26);
	YAPTAZN = y(27);
	REV_mRNA = y(28);
	SRF = y(29);
	CRY = y(30);
	LaminAp = y(31);
	PERCRYnuc_p = y(32);
	YAPTAZnuc = y(33);
	PER_p = y(34);
	mDiaA = y(35);
	LIMKA = y(36);
	mDia = y(37);
	PERCRY = y(38);
	CRY_mRNA = y(39);
	ROCKA = y(40);
	RhoAGTP_MEM = y(41);
	PER = y(42);
	REVnuc = y(43);
	CofilinP = y(44);
	CRY_p = y(45);
	NPC = y(46);
	BMAL1_p = y(47);
	% Constants
	K_MRTF_r46 = p(1);
	Voltage_NM = p(2);
	mlabfix_F_ = p(3);
	CRY_init_uM = p(4);
	YAPTAZN_init_uM = p(5);
	kturnover = p(6);
	PER_init_uM = p(7);
	Fcyto_init_uM = p(8);
	netValence_r59 = p(9);
	netValence_r50 = p(10);
	CRY_p_init_uM = p(11);
	kdrock = p(12);
	PERCRY_init_uM = p(13);
	BMAL1_p_init_uM = p(14);
	ROCKA_init_uM = p(15);
	kCY = p(16);
	netValence_r44 = p(17);
	netValence_r42 = p(18);
	kCN = p(19);
	Positionboolean_init_molecules_um_2 = p(20);
	Gactin_init_uM = p(21);
	ksf = p(22);
	kfc1 = p(23);
	kdmy = p(24);
	krp = p(25);
	PER_p_init_uM = p(26);
	kra = p(27);
	netValence_r22 = p(28);
	netValence_r21 = p(29);
	gamma = p(30);
	SAVol = p(31);
	netValence_r18 = p(32);
	netValence_r17 = p(33);
	netValence_r16 = p(34);
	netValence_r15 = p(35);
	netValence_r14 = p(36);
	Kr_r7 = p(37);
	Kr_r6 = p(38);
	netValence_r12 = p(39);
	Kr_r5 = p(40);
	kcatcof = p(41);
	Kr_r4 = p(42);
	Kr_r3 = p(43);
	PERCRYnuc_init_uM = p(44);
	Kr_r2 = p(45);
	Kr_r1 = p(46);
	Kr_r0 = p(47);
	SAV = p(48);
	K_millivolts_per_volt = p(49);
	mlabfix_K_GHK_ = p(50);
	Size_PM = p(51);
	q_r49 = p(52);
	q_r46 = p(53);
	MyoA_init_uM = p(54);
	kinsolo_MRTF = p(55);
	netValence_r5 = p(56);
	netValence_r4 = p(57);
	netValence_r3 = p(58);
	netValence_r1 = p(59);
	MRTFnuc_init_uM = p(60);
	mlabfix_F_nmol_ = p(61);
	NPCA_init_molecules_um_2 = p(62);
	Size_NM = p(63);
	TEAD_init_uM = p(64);
	q_r27 = p(65);
	kmCof = p(66);
	kmr = p(67);
	kmp = p(68);
	REV_init_uM = p(69);
	mDia_init_uM = p(70);
	kly = p(71);
	MRTFReleaseConst = p(72);
	klr = p(73);
	CofilinNP_init_uM = p(74);
	kll = p(75);
	kout_MRTF = p(76);
	mlabfix_T_ = p(77);
	Size_Cyto = p(78);
	propStiff = p(79);
	Size_Nuc = p(80);
	Clamin = p(81);
	n2 = p(82);
	LaminA_init_molecules_um_2 = p(83);
	n1 = p(84);
	kdep = p(85);
	n_r49 = p(86);
	n_r46 = p(87);
	RhoAGTP_MEM_init_molecules_um_2 = p(88);
	BMAL1_mRNA_init_uM = p(89);
	mlabfix_R_ = p(90);
	kf = p(91);
	CofilinP_init_uM = p(92);
	PERCRY_p_init_uM = p(93);
	LIMK_init_uM = p(94);
	YAPTAZnuc_init_uM = p(95);
	LaminAp_init_molecules_um_2 = p(96);
	REV_mRNA_init_uM = p(97);
	kflaminA = p(98);
	tau = p(99);
	mDiaA_init_uM = p(100);
	BMAL1nuc_p_init_uM = p(101);
	CRY_mRNA_init_uM = p(102);
	BMAL1_init_uM = p(103);
	krNPC = p(104);
	PERCRYnuc_p_init_uM = p(105);
	Fakp_init_uM = p(106);
	kNC = p(107);
	RhoAGDP_init_uM = p(108);
	convertBool = p(109);
	mDiaB = p(110);
	REVnuc_init_uM = p(111);
	Fak_init_uM = p(112);
	Myo_init_uM = p(113);
	kdp = p(114);
	kfkp = p(115);
	kdl = p(116);
	kdf = p(117);
	NPC_init_molecules_um_2 = p(118);
	LIMKA_init_uM = p(119);
	Emol = p(120);
	kfNPC = p(121);
	mlabfix_PI_ = p(122);
	SRF_init_uM = p(123);
	Size_ECM = p(124);
	alpha = p(125);
	ROCKB = p(126);
	PER_mRNA_init_uM = p(127);
	unitconversionfactor = p(128);
	m = p(129);
	kdmdia = p(130);
	h = p(131);
	KMOLE = p(132);
	kin2_MRTF = p(133);
	mlabfix_N_pmol_ = p(134);
	YAPTAZP_init_uM = p(135);
	kinSolo2 = p(136);
	k11 = p(137);
	MRTF_init_uM = p(138);
	C = p(139);
	kout2 = p(140);
	Voltage_PM = p(141);
	K_YTB = p(142);
	epsilon = p(143);
	kin2 = p(144);
	BMAL1nuc_init_uM = p(145);
	ROCK_init_uM = p(146);
	PERCRY_CLOCKBMAL1_init_uM = p(147);
	Kr_r15 = p(148);
	K_MRTF_r49 = p(149);
	Kr_r12 = p(150);
	% Functions
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 .* (KMOLE ^ 1.0));
	nu_dRC = (4.4 .* 0.001 ./ 3600.0);
	K_d_r36 = (0.3 ./ 1000.0);
	K_d_r32 = (0.3 .* 0.001);
	k_sC = (3.2 ./ 3600.0);
	LumpedJ_r44 = (k_sC .* CRY_mRNA .* Size_Cyto ./ KMOLE);
	K_d_r30 = (0.3 .* 0.001);
	k_sP = (1.2 ./ 3600.0);
	LumpedJ_r42 = ((k_sP .* PER_mRNA) .* Size_Cyto ./ KMOLE);
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	k_9 = (0.8 ./ 3600.0);
	k_6 = (0.8 ./ 3600.0);
	k_5 = (0.8 ./ 3600.0);
	k_dnC = (0.02 ./ 3600.0);
	J_r39 = (k_dnC .* CRY);
	K_d_r38 = (0.3 ./ 1000.0);
	k_dn_r38 = (0.02 ./ 3600.0);
	nu_dCC = (1.4 ./ (1000.0 .* 3600.0));
	J_r38 = ((nu_dCC .* CRY_p ./ (K_d_r38 + CRY_p)) + (k_dn_r38 .* CRY_p));
	k_dn_r37 = (0.02 ./ 3600.0);
	J_r37 = (k_dn_r37 .* PERCRY);
	k_2 = (0.4 ./ 3600.0);
	nu_dPCC = (1.4 ./ (1000.0 .* 3600.0));
	k_dn_r36 = (0.02 ./ 3600.0);
	J_r36 = ((nu_dPCC .* PERCRY_p ./ (K_d_r36 + PERCRY_p)) + (k_dn_r36 .* PERCRY_p));
	k_1 = (0.8 ./ 3600.0);
	V_1P = (9.6 ./ (1000.0 .* 3600.0));
	K_dp_r35 = (0.1 ./ 1000.0);
	K_p_r35 = (1.006 ./ 1000.0);
	V_2P = (0.6 ./ (1000.0 .* 3600.0));
	J_r35 = ((V_1P .* PER ./ (K_p_r35 + PER)) - (V_2P .* PER_p ./ (K_dp_r35 + PER_p)));
	V_2C = (0.2 ./ (1000.0 .* 3600.0));
	V_1C = (1.2 ./ (1000.0 .* 3600.0));
	K_p_r34 = (1.006 ./ 1000.0);
	K_dp_r34 = (0.1 ./ 1000.0);
	J_r34 = ((V_1C .* CRY ./ (K_p_r34 + CRY)) - (V_2C .* CRY_p ./ (K_dp_r34 + CRY_p)));
	K_dp_r33 = (0.1 ./ 1000.0);
	K_p_r33 = (1.006 ./ 1000.0);
	V_2PC = (0.2 ./ (1000.0 .* 3600.0));
	V_1PC = (2.4 ./ (1000.0 .* 3600.0));
	J_r33 = ((V_1PC .* PERCRY ./ (K_p_r33 + PERCRY)) - (V_2PC .* PERCRY_p ./ (K_dp_r33 + PERCRY_p)));
	nu_dBN = (3.0 .* 0.001 ./ 3600.0);
	k_dn_r32 = (0.02 ./ 3600.0);
	J_r32 = ((nu_dBN .* BMAL1nuc_p ./ (K_d_r32 + BMAL1nuc_p)) + (k_dn_r32 .* BMAL1nuc_p));
	k_dn_r31 = (0.02 ./ 3600.0);
	J_r31 = (k_dn_r31 .* BMAL1nuc);
	k_dn_r30 = (0.02 ./ 3600.0);
	nu_dBC = (3.0 .* 0.001 ./ 3600.0);
	J_r30 = ((nu_dBC .* BMAL1_p ./ (K_d_r30 + BMAL1_p)) + (k_dn_r30 .* BMAL1_p));
	K_d_r26 = (0.3 .* 0.001);
	K_d_r24 = (0.3 .* 0.001);
	Kr_r9 = kdmy;
	Kf_r9 = (kmr .* (1.0 + (0.5 .* epsilon .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	J_r9 = ((Kf_r9 .* Myo) - (Kr_r9 .* MyoA));
	Kr_r8 = kdl;
	Kf_r8 = (klr .* (1.0 + (0.5 .* tau .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	J_r8 = ((Kf_r8 .* LIMK) - (Kr_r8 .* LIMKA));
	Kf_r7 = kdmdia;
	J_r7 = ((Kf_r7 .* mDiaA) - (Kr_r7 .* mDia));
	Kf_r6 = kdrock;
	J_r6 = ((Kf_r6 .* ROCKA) - (Kr_r6 .* ROCK));
	Kf_r5 = (kmp .* RhoAGTP_MEM);
	J_r5 = ((Kf_r5 .* mDia) - (Kr_r5 .* mDiaA));
	Kf_r4 = (krp .* RhoAGTP_MEM);
	J_r4 = ((Kf_r4 .* ROCK) - (Kr_r4 .* ROCKA));
	Kf_r3 = kdp;
	J_r3 = ((Kf_r3 .* RhoAGTP_MEM) - (Kr_r3 .* RhoAGDP));
	k_dn_r29 = (0.02 ./ 3600.0);
	J_r29 = (k_dn_r29 .* BMAL1);
	k_dmB = (0.02 ./ 3600.0);
	nu_mB = (1.3 .* 0.001 ./ 3600.0);
	K_mB = (0.4 .* 0.001);
	J_r28 = (((nu_mB .* BMAL1_mRNA) ./ (K_mB + BMAL1_mRNA)) + (k_dmB .* BMAL1_mRNA));
	Kf_r2 = kf;
	J_r2 = ((Kf_r2 .* Fak) - (Kr_r2 .* Fakp));
	nu_sB = (1.8 .* 0.001 ./ 3600.0);
	nu_YTB = (0.3 .* 0.001 ./ 3600.0);
	K_IB = (2.2 .* 0.001);
	J_r27 = (((nu_sB .* (K_IB ^ m)) ./ ((K_IB ^ m) + (REVnuc ^ m))) + ((nu_YTB .* (YAPTAZnuc ^ q_r27)) ./ ((K_YTB ^ q_r27) + (YAPTAZnuc ^ q_r27))));
	Kf_r1 = (kfkp .* (1.0 + (gamma .* (Fakp ^ n2))) .* unitconversionfactor .* SAV);
	J_r1 = ((Kf_r1 .* RhoAGDP) - (Kr_r1 .* RhoAGTP_MEM));
	k_dn_r26 = (0.02 ./ 3600.0);
	nu_dRN = (0.8 .* 0.001 ./ 3600.0);
	J_r26 = (((nu_dRN .* REVnuc) ./ (K_d_r26 + REVnuc)) - (k_dn_r26 .* REVnuc));
	Kf_r0 = kdf;
	J_r0 = ((Kf_r0 .* Fakp) - (Kr_r0 .* Fak));
	k_dmr = (0.02 ./ 3600.0);
	nu_mR = (1.6 .* 0.001 ./ 3600.0);
	K_mR = (0.4 .* 0.001);
	J_r25 = (((nu_mR .* REV_mRNA) ./ (K_mR + REV_mRNA)) + (k_dmr .* REV_mRNA));
	k_dn_r24 = (0.02 ./ 3600.0);
	J_r24 = (((nu_dRC .* REV) ./ (K_d_r24 + REV)) + (k_dn_r24 .* REV));
	nu_sR = (1.6 .* 0.001 ./ 3600.0);
	K_AR = (0.6 .* 0.001);
	J_r23 = ((nu_sR .* (BMAL1nuc ^ h)) ./ ((K_AR ^ h) + (BMAL1nuc ^ h)));
	V_2B = (0.2 .* 0.001 ./ 3600.0);
	V_1B = (1.4 .* 0.001 ./ 3600.0);
	K_dp_r20 = (0.1 .* 0.001);
	K_p_r20 = (1.006 .* 0.001);
	J_r20 = ((V_1B .* BMAL1 ./ (K_p_r20 + BMAL1)) - (V_2B .* BMAL1_p ./ (K_dp_r20 + BMAL1_p)));
	nu_dPC = (3.4 ./ (1000.0 .* 3600.0));
	k_10 = (0.4 ./ 3600.0);
	Kr_r22 = k_10;
	Kf_r22 = k_9;
	LumpedJ_r22 = (((Kf_r22 .* REV) - (Kr_r22 .* REVnuc)) .* Size_Cyto ./ KMOLE);
	LumpedJ_r21 = (((k_5 .* BMAL1) - (k_6 .* BMAL1nuc)) .* Size_Cyto ./ KMOLE);
	K_p_r52 = (1.006 ./ 1000.0);
	V_3B = (1.4 .* 0.001 ./ 3600.0);
	K_p_r19 = (1.006 .* 0.001);
	K_dp_r19 = (0.1 .* 0.001);
	V_4B = (0.4 .* 0.001 ./ 3600.0);
	J_r19 = ((V_3B .* BMAL1nuc ./ (K_p_r19 + BMAL1nuc)) - (V_4B .* BMAL1nuc_p ./ (K_dp_r19 + BMAL1nuc_p)));
	KFlux_NM_Nuc = (Size_NM ./ Size_Nuc);
	Kr_r16 = kout2;
	Kf_r16 = ((kin2 .* NPCA) + kinSolo2);
	J_r16 = ((Kf_r16 .* YAPTAZN) - (Kr_r16 .* YAPTAZnuc));
	Ecytosol = (propStiff .* (Fcyto ^ n1));
	Kf_r15 = (kflaminA .* Ecytosol ./ (Clamin + Ecytosol));
	J_r15 = ((Kf_r15 .* LaminAp) - (Kr_r15 .* LaminA));
	Kr_r14 = krNPC;
	Kf_r14 = (kfNPC .* MyoA .* Fcyto .* LaminA);
	J_r14 = ((Kf_r14 .* NPC) - (Kr_r14 .* NPCA));
	Kr_r13 = kNC;
	Kf_r13 = (kCN + (kCY .* MyoA .* Fcyto));
	J_r13 = ((Kf_r13 .* YAPTAZP) - (Kr_r13 .* YAPTAZN));
	Kf_r12 = (ksf .* Emol .* unitconversionfactor .* SAVol .* Positionboolean .* convertBool ./ (C + Emol));
	J_r12 = ((Kf_r12 .* Fak) - (Kr_r12 .* Fakp));
	Kr_r11 = (kdep + (kfc1 .* CofilinNP));
	Kf_r11 = (kra .* (1.0 + (0.5 .* alpha .* (1.0 + tanh((20.0 .* (mDiaA - mDiaB)))) .* mDiaA)));
	J_r11 = ((Kf_r11 .* Gactin) - (Kr_r11 .* Fcyto));
	Kr_r10 = (kcatcof .* LIMKA ./ (kmCof + CofilinNP));
	Kf_r10 = kturnover;
	J_r10 = ((Kf_r10 .* CofilinP) - (Kr_r10 .* CofilinNP));
	k_sB = (0.32 ./ 3600.0);
	LumpedJ_r18 = (k_sB .* BMAL1_mRNA .* Size_Cyto ./ KMOLE);
	k_sR = (1.7 ./ 3600.0);
	LumpedJ_r17 = (k_sR .* REV_mRNA .* Size_Cyto ./ KMOLE);
	K_AP = (0.6 ./ 1000.0);
	K_AC = (0.6 ./ 1000.0);
	nu_MRTF_r49 = (1.0 ./ (1000.0 .* 3600.0));
	nu_MRTF_r46 = (1.0 ./ (1000.0 .* 3600.0));
	k_4_r51 = (0.4 ./ 3600.0);
	nu_sP = (2.4 ./ (1000.0 .* 3600.0));
	k_4_r45 = (0.4 ./ 3600.0);
	k_4_r43 = (0.4 ./ 3600.0);
	nu_sC = (2.2 ./ (1000.0 .* 3600.0));
	k_3_r51 = (0.8 .* 1000.0 ./ 3600.0);
	nu_dPCN = (1.4 ./ (1000.0 .* 3600.0));
	k_3_r45 = (0.8 .* 1000.0 ./ 3600.0);
	k_3_r43 = (0.8 .* 1000.0 ./ 3600.0);
	nu_dIN = (1.6 ./ (1000.0 .* 3600.0));
	K_mP = (0.3 ./ 1000.0);
	K_mC = (0.4 ./ 1000.0);
	nu_mP = (2.2 ./ (1000.0 .* 3600.0));
	nu_mC = (2.0 ./ (1000.0 .* 3600.0));
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	k_dn_r58 = (0.02 ./ 3600.0);
	k_dn_r54 = (0.02 ./ 3600.0);
	k_dn_r53 = (0.02 ./ 3600.0);
	k_dn_r41 = (0.02 ./ 3600.0);
	k_dn_r40 = (0.02 ./ 3600.0);
	K_dp_r52 = (0.1 ./ 1000.0);
	Kf_r59 = ((1.0 ./ (1.0 + ((Gactin ./ MRTFReleaseConst) ^ 2.0))) .* (kinsolo_MRTF + (kin2_MRTF .* NPCA)));
	V_4PC = (0.2 ./ (1000.0 .* 3600.0));
	Kr_r59 = kout_MRTF;
	k_dmp = (0.02 ./ 3600.0);
	k_dmc = (0.02 ./ 3600.0);
	k_8_r57 = (0.2 ./ 3600.0);
	K_d_r58 = (0.3 ./ 1000.0);
	k_8_r56 = (0.2 ./ 3600.0);
	k_8_r55 = (0.2 ./ 3600.0);
	K_d_r54 = (0.3 ./ 1000.0);
	J_r59 = ((Kf_r59 .* MRTF) - (Kr_r59 .* MRTFnuc));
	V_3PC = (2.4 ./ (1000.0 .* 3600.0));
	J_r58 = ((nu_dIN .* PERCRY_CLOCKBMAL1 ./ (K_d_r58 + PERCRY_CLOCKBMAL1)) + (k_dn_r58 .* PERCRY_CLOCKBMAL1));
	k_7_r57 = (1.0 .* 1000.0 ./ 3600.0);
	J_r57 = ( - (k_8_r57 .* PERCRY_CLOCKBMAL1) + (k_7_r57 .* BMAL1nuc .* PERCRYnuc));
	k_7_r56 = (1.0 .* 1000.0 ./ 3600.0);
	J_r56 = ((k_8_r56 .* PERCRY_CLOCKBMAL1) - (k_7_r56 .* BMAL1nuc .* PERCRYnuc));
	k_7_r55 = (1.0 .* 1000.0 ./ 3600.0);
	J_r55 = ((k_8_r55 .* PERCRY_CLOCKBMAL1) - (k_7_r55 .* BMAL1nuc .* PERCRYnuc));
	J_r54 = ((nu_dPCN .* PERCRYnuc_p ./ (K_d_r54 + PERCRYnuc_p)) + (k_dn_r54 .* PERCRYnuc_p));
	J_r53 = (k_dn_r53 .* PERCRYnuc);
	J_r52 = ((V_3PC .* PERCRYnuc ./ (K_p_r52 + PERCRYnuc)) - (V_4PC .* PERCRYnuc_p ./ (K_dp_r52 + PERCRYnuc_p)));
	J_r51 = ( - (k_4_r51 .* PERCRY) + (k_3_r51 .* PER .* CRY));
	K_d_r40 = (0.3 ./ 1000.0);
	LumpedJ_r50 = (((k_1 .* PERCRY) - (k_2 .* PERCRYnuc)) .* Size_Cyto ./ KMOLE);
	J_r49 = ((nu_sC .* (BMAL1nuc ^ n_r49) ./ ((K_AC ^ n_r49) + (BMAL1nuc ^ n_r49))) + (nu_MRTF_r49 .* (MRTFnuc ^ q_r49) ./ ((K_MRTF_r49 ^ q_r49) + (MRTFnuc ^ q_r49))));
	J_r48 = ((nu_mC .* CRY_mRNA ./ (K_mC + CRY_mRNA)) + (k_dmc .* CRY_mRNA));
	J_r47 = ((nu_mP .* PER_mRNA ./ (K_mP + PER_mRNA)) + (k_dmp .* PER_mRNA));
	J_r46 = ((nu_sP .* (BMAL1nuc ^ n_r46) ./ ((K_AP ^ n_r46) + (BMAL1nuc ^ n_r46))) + (nu_MRTF_r46 .* (MRTFnuc ^ q_r46) ./ ((K_MRTF_r46 ^ q_r46) + (MRTFnuc ^ q_r46))));
	J_r45 = ((k_4_r45 .* PERCRY) - (k_3_r45 .* PER .* CRY));
	J_r43 = ((k_4_r43 .* PERCRY) - (k_3_r43 .* PER .* CRY));
	J_r41 = (k_dn_r41 .* PER);
	J_r40 = ((nu_dPC .* PER_p ./ (K_d_r40 + PER_p)) + (k_dn_r40 .* PER_p));
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

	rowValue = [MRTF PERCRYnuc LaminA MyoA BMAL1 Fakp TEAD PERCRY_CLOCKBMAL1 Gactin Fak ROCK BMAL1_mRNA CofilinNP REV BMAL1nuc NPCA PER_mRNA Positionboolean LIMK RhoAGDP MRTFnuc Fcyto BMAL1nuc_p Myo PERCRY_p YAPTAZP YAPTAZN REV_mRNA SRF CRY LaminAp PERCRYnuc_p YAPTAZnuc PER_p mDiaA LIMKA mDia PERCRY CRY_mRNA ROCKA RhoAGTP_MEM PER REVnuc CofilinP CRY_p NPC BMAL1_p UnitFactor_uM_um3_molecules_neg_1 nu_dRC K_d_r36 K_d_r32 k_sC LumpedJ_r44 K_d_r30 k_sP LumpedJ_r42 KFlux_PM_Cyto k_9 k_6 k_5 k_dnC J_r39 K_d_r38 k_dn_r38 nu_dCC J_r38 k_dn_r37 J_r37 k_2 nu_dPCC k_dn_r36 J_r36 k_1 V_1P K_dp_r35 K_p_r35 V_2P J_r35 V_2C V_1C K_p_r34 K_dp_r34 J_r34 K_dp_r33 K_p_r33 V_2PC V_1PC J_r33 nu_dBN k_dn_r32 J_r32 k_dn_r31 J_r31 k_dn_r30 nu_dBC J_r30 K_d_r26 K_d_r24 Kr_r9 Kf_r9 J_r9 Kr_r8 Kf_r8 J_r8 Kf_r7 J_r7 Kf_r6 J_r6 Kf_r5 J_r5 Kf_r4 J_r4 Kf_r3 J_r3 k_dn_r29 J_r29 k_dmB nu_mB K_mB J_r28 Kf_r2 J_r2 nu_sB nu_YTB K_IB J_r27 Kf_r1 J_r1 k_dn_r26 nu_dRN J_r26 Kf_r0 J_r0 k_dmr nu_mR K_mR J_r25 k_dn_r24 J_r24 nu_sR K_AR J_r23 V_2B V_1B K_dp_r20 K_p_r20 J_r20 nu_dPC k_10 Kr_r22 Kf_r22 LumpedJ_r22 LumpedJ_r21 K_p_r52 V_3B K_p_r19 K_dp_r19 V_4B J_r19 KFlux_NM_Nuc Kr_r16 Kf_r16 J_r16 Ecytosol Kf_r15 J_r15 Kr_r14 Kf_r14 J_r14 Kr_r13 Kf_r13 J_r13 Kf_r12 J_r12 Kr_r11 Kf_r11 J_r11 Kr_r10 Kf_r10 J_r10 k_sB LumpedJ_r18 k_sR LumpedJ_r17 K_AP K_AC nu_MRTF_r49 nu_MRTF_r46 k_4_r51 nu_sP k_4_r45 k_4_r43 nu_sC k_3_r51 nu_dPCN k_3_r45 k_3_r43 nu_dIN K_mP K_mC nu_mP nu_mC KFlux_NM_Cyto k_dn_r58 k_dn_r54 k_dn_r53 k_dn_r41 k_dn_r40 K_dp_r52 Kf_r59 V_4PC Kr_r59 k_dmp k_dmc k_8_r57 K_d_r58 k_8_r56 k_8_r55 K_d_r54 J_r59 V_3PC J_r58 k_7_r57 J_r57 k_7_r56 J_r56 k_7_r55 J_r55 J_r54 J_r53 J_r52 J_r51 K_d_r40 LumpedJ_r50 J_r49 J_r48 J_r47 J_r46 J_r45 J_r43 J_r41 J_r40 YAPTAZratio FakConservation RhoAConservation LIMKConservation YAPTAZConservation MyoConservation mDiaConservation ROCKConservation CofilinConservation LaminAConservation NPCConservation ActinConservation];
end

% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0,couplingParam)
	% State Variables
	MRTF = y(1);
	PERCRYnuc = y(2);
	LaminA = y(3);
	MyoA = y(4);
	BMAL1 = y(5);
	Fakp = y(6);
	TEAD = y(7);
	PERCRY_CLOCKBMAL1 = y(8);
	Gactin = y(9);
	Fak = y(10);
	ROCK = y(11);
	BMAL1_mRNA = y(12);
	CofilinNP = y(13);
	REV = y(14);
	BMAL1nuc = y(15);
	NPCA = y(16);
	PER_mRNA = y(17);
	Positionboolean = y(18);
	LIMK = y(19);
	RhoAGDP = y(20);
	MRTFnuc = y(21);
	Fcyto = y(22);
	BMAL1nuc_p = y(23);
	Myo = y(24);
	PERCRY_p = y(25);
	YAPTAZP = y(26);
	YAPTAZN = y(27);
	REV_mRNA = y(28);
	SRF = y(29);
	CRY = y(30);
	LaminAp = y(31);
	PERCRYnuc_p = y(32);
	YAPTAZnuc = y(33);
	PER_p = y(34);
	mDiaA = y(35);
	LIMKA = y(36);
	mDia = y(37);
	PERCRY = y(38);
	CRY_mRNA = y(39);
	ROCKA = y(40);
	RhoAGTP_MEM = y(41);
	PER = y(42);
	REVnuc = y(43);
	CofilinP = y(44);
	CRY_p = y(45);
	NPC = y(46);
	BMAL1_p = y(47);
	% Constants
	K_MRTF_r46 = p(1);
	Voltage_NM = p(2);
	mlabfix_F_ = p(3);
	CRY_init_uM = p(4);
	YAPTAZN_init_uM = p(5);
	kturnover = p(6);
	PER_init_uM = p(7);
	Fcyto_init_uM = p(8);
	netValence_r59 = p(9);
	netValence_r50 = p(10);
	CRY_p_init_uM = p(11);
	kdrock = p(12);
	PERCRY_init_uM = p(13);
	BMAL1_p_init_uM = p(14);
	ROCKA_init_uM = p(15);
	kCY = p(16);
	netValence_r44 = p(17);
	netValence_r42 = p(18);
	kCN = p(19);
	Positionboolean_init_molecules_um_2 = p(20);
	Gactin_init_uM = p(21);
	ksf = p(22);
	kfc1 = p(23);
	kdmy = p(24);
	krp = p(25);
	PER_p_init_uM = p(26);
	kra = p(27);
	netValence_r22 = p(28);
	netValence_r21 = p(29);
	gamma = p(30);
	SAVol = p(31);
	netValence_r18 = p(32);
	netValence_r17 = p(33);
	netValence_r16 = p(34);
	netValence_r15 = p(35);
	netValence_r14 = p(36);
	Kr_r7 = p(37);
	Kr_r6 = p(38);
	netValence_r12 = p(39);
	Kr_r5 = p(40);
	kcatcof = p(41);
	Kr_r4 = p(42);
	Kr_r3 = p(43);
	PERCRYnuc_init_uM = p(44);
	Kr_r2 = p(45);
	Kr_r1 = p(46);
	Kr_r0 = p(47);
	SAV = p(48);
	K_millivolts_per_volt = p(49);
	mlabfix_K_GHK_ = p(50);
	Size_PM = p(51);
	q_r49 = p(52);
	q_r46 = p(53);
	MyoA_init_uM = p(54);
	kinsolo_MRTF = p(55);
	netValence_r5 = p(56);
	netValence_r4 = p(57);
	netValence_r3 = p(58);
	netValence_r1 = p(59);
	MRTFnuc_init_uM = p(60);
	mlabfix_F_nmol_ = p(61);
	NPCA_init_molecules_um_2 = p(62);
	Size_NM = p(63);
	TEAD_init_uM = p(64);
	q_r27 = p(65);
	kmCof = p(66);
	kmr = p(67);
	kmp = p(68);
	REV_init_uM = p(69);
	mDia_init_uM = p(70);
	kly = p(71);
	MRTFReleaseConst = p(72);
	klr = p(73);
	CofilinNP_init_uM = p(74);
	kll = p(75);
	kout_MRTF = p(76);
	mlabfix_T_ = p(77);
	Size_Cyto = p(78);
	propStiff = p(79);
	Size_Nuc = p(80);
	Clamin = p(81);
	n2 = p(82);
	LaminA_init_molecules_um_2 = p(83);
	n1 = p(84);
	kdep = p(85);
	n_r49 = p(86);
	n_r46 = p(87);
	RhoAGTP_MEM_init_molecules_um_2 = p(88);
	BMAL1_mRNA_init_uM = p(89);
	mlabfix_R_ = p(90);
	kf = p(91);
	CofilinP_init_uM = p(92);
	PERCRY_p_init_uM = p(93);
	LIMK_init_uM = p(94);
	YAPTAZnuc_init_uM = p(95);
	LaminAp_init_molecules_um_2 = p(96);
	REV_mRNA_init_uM = p(97);
	kflaminA = p(98);
	tau = p(99);
	mDiaA_init_uM = p(100);
	BMAL1nuc_p_init_uM = p(101);
	CRY_mRNA_init_uM = p(102);
	BMAL1_init_uM = p(103);
	krNPC = p(104);
	PERCRYnuc_p_init_uM = p(105);
	Fakp_init_uM = p(106);
	kNC = p(107);
	RhoAGDP_init_uM = p(108);
	convertBool = p(109);
	mDiaB = p(110);
	REVnuc_init_uM = p(111);
	Fak_init_uM = p(112);
	Myo_init_uM = p(113);
	kdp = p(114);
	kfkp = p(115);
	kdl = p(116);
	kdf = p(117);
	NPC_init_molecules_um_2 = p(118);
	LIMKA_init_uM = p(119);
	Emol = p(120);
	kfNPC = p(121);
	mlabfix_PI_ = p(122);
	SRF_init_uM = p(123);
	Size_ECM = p(124);
	alpha = p(125);
	ROCKB = p(126);
	PER_mRNA_init_uM = p(127);
	unitconversionfactor = p(128);
	m = p(129);
	kdmdia = p(130);
	h = p(131);
	KMOLE = p(132);
	kin2_MRTF = p(133);
	mlabfix_N_pmol_ = p(134);
	YAPTAZP_init_uM = p(135);
	kinSolo2 = p(136);
	k11 = p(137);
	MRTF_init_uM = p(138);
	C = p(139);
	kout2 = p(140);
	Voltage_PM = p(141);
	K_YTB = p(142);
	epsilon = p(143);
	kin2 = p(144);
	BMAL1nuc_init_uM = p(145);
	ROCK_init_uM = p(146);
	PERCRY_CLOCKBMAL1_init_uM = p(147);
	Kr_r15 = p(148);
	K_MRTF_r49 = p(149);
	Kr_r12 = p(150);
	% Functions
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 .* (KMOLE ^ 1.0));
	nu_dRC = (4.4 .* 0.001 ./ 3600.0);
	K_d_r36 = (0.3 ./ 1000.0);
	K_d_r32 = (0.3 .* 0.001);
	k_sC = (3.2 ./ 3600.0);
	LumpedJ_r44 = (k_sC .* CRY_mRNA .* Size_Cyto ./ KMOLE);
	K_d_r30 = (0.3 .* 0.001);
	k_sP = (1.2 ./ 3600.0);
	LumpedJ_r42 = ((k_sP .* PER_mRNA) .* Size_Cyto ./ KMOLE);
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	k_9 = (0.8 ./ 3600.0);
	k_6 = (0.8 ./ 3600.0);
	k_5 = (0.8 ./ 3600.0);
	k_dnC = (0.02 ./ 3600.0);
	J_r39 = (k_dnC .* CRY);
	K_d_r38 = (0.3 ./ 1000.0);
	k_dn_r38 = (0.02 ./ 3600.0);
	nu_dCC = (1.4 ./ (1000.0 .* 3600.0));
	J_r38 = ((nu_dCC .* CRY_p ./ (K_d_r38 + CRY_p)) + (k_dn_r38 .* CRY_p));
	k_dn_r37 = (0.02 ./ 3600.0);
	J_r37 = (k_dn_r37 .* PERCRY);
	k_2 = (0.4 ./ 3600.0);
	nu_dPCC = (1.4 ./ (1000.0 .* 3600.0));
	k_dn_r36 = (0.02 ./ 3600.0);
	J_r36 = ((nu_dPCC .* PERCRY_p ./ (K_d_r36 + PERCRY_p)) + (k_dn_r36 .* PERCRY_p));
	k_1 = (0.8 ./ 3600.0);
	V_1P = (9.6 ./ (1000.0 .* 3600.0));
	K_dp_r35 = (0.1 ./ 1000.0);
	K_p_r35 = (1.006 ./ 1000.0);
	V_2P = (0.6 ./ (1000.0 .* 3600.0));
	J_r35 = ((V_1P .* PER ./ (K_p_r35 + PER)) - (V_2P .* PER_p ./ (K_dp_r35 + PER_p)));
	V_2C = (0.2 ./ (1000.0 .* 3600.0));
	V_1C = (1.2 ./ (1000.0 .* 3600.0));
	K_p_r34 = (1.006 ./ 1000.0);
	K_dp_r34 = (0.1 ./ 1000.0);
	J_r34 = ((V_1C .* CRY ./ (K_p_r34 + CRY)) - (V_2C .* CRY_p ./ (K_dp_r34 + CRY_p)));
	K_dp_r33 = (0.1 ./ 1000.0);
	K_p_r33 = (1.006 ./ 1000.0);
	V_2PC = (0.2 ./ (1000.0 .* 3600.0));
	V_1PC = (2.4 ./ (1000.0 .* 3600.0));
	J_r33 = ((V_1PC .* PERCRY ./ (K_p_r33 + PERCRY)) - (V_2PC .* PERCRY_p ./ (K_dp_r33 + PERCRY_p)));
	nu_dBN = (3.0 .* 0.001 ./ 3600.0);
	k_dn_r32 = (0.02 ./ 3600.0);
	J_r32 = ((nu_dBN .* BMAL1nuc_p ./ (K_d_r32 + BMAL1nuc_p)) + (k_dn_r32 .* BMAL1nuc_p));
	k_dn_r31 = (0.02 ./ 3600.0);
	J_r31 = (k_dn_r31 .* BMAL1nuc);
	k_dn_r30 = (0.02 ./ 3600.0);
	nu_dBC = (3.0 .* 0.001 ./ 3600.0);
	J_r30 = ((nu_dBC .* BMAL1_p ./ (K_d_r30 + BMAL1_p)) + (k_dn_r30 .* BMAL1_p));
	K_d_r26 = (0.3 .* 0.001);
	K_d_r24 = (0.3 .* 0.001);
	Kr_r9 = kdmy;
	Kf_r9 = (kmr .* (1.0 + (0.5 .* epsilon .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	J_r9 = ((Kf_r9 .* Myo) - (Kr_r9 .* MyoA));
	Kr_r8 = kdl;
	Kf_r8 = (klr .* (1.0 + (0.5 .* tau .* (1.0 + tanh((20.0 .* (ROCKA - ROCKB)))) .* ROCKA)));
	J_r8 = ((Kf_r8 .* LIMK) - (Kr_r8 .* LIMKA));
	Kf_r7 = kdmdia;
	J_r7 = ((Kf_r7 .* mDiaA) - (Kr_r7 .* mDia));
	Kf_r6 = kdrock;
	J_r6 = ((Kf_r6 .* ROCKA) - (Kr_r6 .* ROCK));
	Kf_r5 = (kmp .* RhoAGTP_MEM);
	J_r5 = ((Kf_r5 .* mDia) - (Kr_r5 .* mDiaA));
	Kf_r4 = (krp .* RhoAGTP_MEM);
	J_r4 = ((Kf_r4 .* ROCK) - (Kr_r4 .* ROCKA));
	Kf_r3 = kdp;
	J_r3 = ((Kf_r3 .* RhoAGTP_MEM) - (Kr_r3 .* RhoAGDP));
	k_dn_r29 = (0.02 ./ 3600.0);
	J_r29 = (k_dn_r29 .* BMAL1);
	k_dmB = (0.02 ./ 3600.0);
	nu_mB = (1.3 .* 0.001 ./ 3600.0);
	K_mB = (0.4 .* 0.001);
	J_r28 = (((nu_mB .* BMAL1_mRNA) ./ (K_mB + BMAL1_mRNA)) + (k_dmB .* BMAL1_mRNA));
	Kf_r2 = kf;
	J_r2 = ((Kf_r2 .* Fak) - (Kr_r2 .* Fakp));
	nu_sB = (1.8 .* 0.001 ./ 3600.0);
	nu_YTB = couplingParam(1,1);%(0.3 .* 0.001 ./ 3600.0);
    nu_YTP = couplingParam(1,2);
    nu_YTC = couplingParam(1,3);
    nu_YTR = couplingParam(1,4);
    nu_MRTF_B = couplingParam(2,1);%(1.0 ./ (1000.0 .* 3600.0));
	nu_MRTF_P = couplingParam(2,2);%(1.0 ./ (1000.0 .* 3600.0));
    nu_MRTF_C = couplingParam(2,3);
    nu_MRTF_R = couplingParam(2,4);

	K_IB = (2.2 .* 0.001);
    YAPTAZHillFcn = ((YAPTAZnuc ^ q_r27)) ./ ((K_YTB ^ q_r27) + (YAPTAZnuc ^ q_r27));
    MRTFHillFcn = (MRTFnuc ^ q_r49) ./ ((K_MRTF_r49 ^ q_r49) + (MRTFnuc ^ q_r49));
	J_r27 = (((nu_sB .* (K_IB ^ m)) ./ ((K_IB ^ m) + (REVnuc ^ m))) + (nu_YTB .* YAPTAZHillFcn) + nu_MRTF_B.*MRTFHillFcn);
	Kf_r1 = (kfkp .* (1.0 + (gamma .* (Fakp ^ n2))) .* unitconversionfactor .* SAV);
	J_r1 = ((Kf_r1 .* RhoAGDP) - (Kr_r1 .* RhoAGTP_MEM));
	k_dn_r26 = (0.02 ./ 3600.0);
	nu_dRN = (0.8 .* 0.001 ./ 3600.0);
	J_r26 = (((nu_dRN .* REVnuc) ./ (K_d_r26 + REVnuc)) - (k_dn_r26 .* REVnuc));
	Kf_r0 = kdf;
	J_r0 = ((Kf_r0 .* Fakp) - (Kr_r0 .* Fak));
	k_dmr = (0.02 ./ 3600.0);
	nu_mR = (1.6 .* 0.001 ./ 3600.0);
	K_mR = (0.4 .* 0.001);
	J_r25 = (((nu_mR .* REV_mRNA) ./ (K_mR + REV_mRNA)) + (k_dmr .* REV_mRNA)) + nu_YTR.*YAPTAZHillFcn + nu_MRTF_R.*MRTFHillFcn;
	k_dn_r24 = (0.02 ./ 3600.0);
	J_r24 = (((nu_dRC .* REV) ./ (K_d_r24 + REV)) + (k_dn_r24 .* REV));
	nu_sR = (1.6 .* 0.001 ./ 3600.0);
	K_AR = (0.6 .* 0.001);
	J_r23 = ((nu_sR .* (BMAL1nuc ^ h)) ./ ((K_AR ^ h) + (BMAL1nuc ^ h)));
	V_2B = (0.2 .* 0.001 ./ 3600.0);
	V_1B = (1.4 .* 0.001 ./ 3600.0);
	K_dp_r20 = (0.1 .* 0.001);
	K_p_r20 = (1.006 .* 0.001);
	J_r20 = ((V_1B .* BMAL1 ./ (K_p_r20 + BMAL1)) - (V_2B .* BMAL1_p ./ (K_dp_r20 + BMAL1_p)));
	nu_dPC = (3.4 ./ (1000.0 .* 3600.0));
	k_10 = (0.4 ./ 3600.0);
	Kr_r22 = k_10;
	Kf_r22 = k_9;
	LumpedJ_r22 = (((Kf_r22 .* REV) - (Kr_r22 .* REVnuc)) .* Size_Cyto ./ KMOLE);
	LumpedJ_r21 = (((k_5 .* BMAL1) - (k_6 .* BMAL1nuc)) .* Size_Cyto ./ KMOLE);
	K_p_r52 = (1.006 ./ 1000.0);
	V_3B = (1.4 .* 0.001 ./ 3600.0);
	K_p_r19 = (1.006 .* 0.001);
	K_dp_r19 = (0.1 .* 0.001);
	V_4B = (0.4 .* 0.001 ./ 3600.0);
	J_r19 = ((V_3B .* BMAL1nuc ./ (K_p_r19 + BMAL1nuc)) - (V_4B .* BMAL1nuc_p ./ (K_dp_r19 + BMAL1nuc_p)));
	KFlux_NM_Nuc = (Size_NM ./ Size_Nuc);
	Kr_r16 = kout2;
	Kf_r16 = ((kin2 .* NPCA) + kinSolo2);
	J_r16 = ((Kf_r16 .* YAPTAZN) - (Kr_r16 .* YAPTAZnuc));
	Ecytosol = (propStiff .* (Fcyto ^ n1));
	Kf_r15 = (kflaminA .* Ecytosol ./ (Clamin + Ecytosol));
	J_r15 = ((Kf_r15 .* LaminAp) - (Kr_r15 .* LaminA));
	Kr_r14 = krNPC;
	Kf_r14 = (kfNPC .* MyoA .* Fcyto .* LaminA);
	J_r14 = ((Kf_r14 .* NPC) - (Kr_r14 .* NPCA));
	Kr_r13 = kNC;
	Kf_r13 = (kCN + (kCY .* MyoA .* Fcyto));
	J_r13 = ((Kf_r13 .* YAPTAZP) - (Kr_r13 .* YAPTAZN));
	Kf_r12 = (ksf .* Emol .* unitconversionfactor .* SAVol .* Positionboolean .* convertBool ./ (C + Emol));
	J_r12 = ((Kf_r12 .* Fak) - (Kr_r12 .* Fakp));
	Kr_r11 = (kdep + (kfc1 .* CofilinNP));
	Kf_r11 = (kra .* (1.0 + (0.5 .* alpha .* (1.0 + tanh((20.0 .* (mDiaA - mDiaB)))) .* mDiaA)));
	J_r11 = ((Kf_r11 .* Gactin) - (Kr_r11 .* Fcyto));
	Kr_r10 = (kcatcof .* LIMKA ./ (kmCof + CofilinNP));
	Kf_r10 = kturnover;
	J_r10 = ((Kf_r10 .* CofilinP) - (Kr_r10 .* CofilinNP));
	k_sB = (0.32 ./ 3600.0);
	LumpedJ_r18 = (k_sB .* BMAL1_mRNA .* Size_Cyto ./ KMOLE);
	k_sR = (1.7 ./ 3600.0);
	LumpedJ_r17 = (k_sR .* REV_mRNA .* Size_Cyto ./ KMOLE);
	K_AP = (0.6 ./ 1000.0);
	K_AC = (0.6 ./ 1000.0);
	
	k_4_r51 = (0.4 ./ 3600.0);
	nu_sP = (2.4 ./ (1000.0 .* 3600.0));
	k_4_r45 = (0.4 ./ 3600.0);
	k_4_r43 = (0.4 ./ 3600.0);
	nu_sC = (2.2 ./ (1000.0 .* 3600.0));
	k_3_r51 = (0.8 .* 1000.0 ./ 3600.0);
	nu_dPCN = (1.4 ./ (1000.0 .* 3600.0));
	k_3_r45 = (0.8 .* 1000.0 ./ 3600.0);
	k_3_r43 = (0.8 .* 1000.0 ./ 3600.0);
	nu_dIN = (1.6 ./ (1000.0 .* 3600.0));
	K_mP = (0.3 ./ 1000.0);
	K_mC = (0.4 ./ 1000.0);
	nu_mP = (2.2 ./ (1000.0 .* 3600.0));
	nu_mC = (2.0 ./ (1000.0 .* 3600.0));
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	k_dn_r58 = (0.02 ./ 3600.0);
	k_dn_r54 = (0.02 ./ 3600.0);
	k_dn_r53 = (0.02 ./ 3600.0);
	k_dn_r41 = (0.02 ./ 3600.0);
	k_dn_r40 = (0.02 ./ 3600.0);
	K_dp_r52 = (0.1 ./ 1000.0);
	Kf_r59 = ((1.0 ./ (1.0 + ((Gactin ./ MRTFReleaseConst) ^ 2.0))) .* (kinsolo_MRTF + (kin2_MRTF .* NPCA)));
	V_4PC = (0.2 ./ (1000.0 .* 3600.0));
	Kr_r59 = kout_MRTF;
	k_dmp = (0.02 ./ 3600.0);
	k_dmc = (0.02 ./ 3600.0);
	k_8_r57 = (0.2 ./ 3600.0);
	K_d_r58 = (0.3 ./ 1000.0);
	k_8_r56 = (0.2 ./ 3600.0);
	k_8_r55 = (0.2 ./ 3600.0);
	K_d_r54 = (0.3 ./ 1000.0);
	J_r59 = ((Kf_r59 .* MRTF) - (Kr_r59 .* MRTFnuc));
	V_3PC = (2.4 ./ (1000.0 .* 3600.0));
	J_r58 = ((nu_dIN .* PERCRY_CLOCKBMAL1 ./ (K_d_r58 + PERCRY_CLOCKBMAL1)) + (k_dn_r58 .* PERCRY_CLOCKBMAL1));
	k_7_r57 = (1.0 .* 1000.0 ./ 3600.0);
	J_r57 = ( - (k_8_r57 .* PERCRY_CLOCKBMAL1) + (k_7_r57 .* BMAL1nuc .* PERCRYnuc));
	k_7_r56 = (1.0 .* 1000.0 ./ 3600.0);
	J_r56 = ((k_8_r56 .* PERCRY_CLOCKBMAL1) - (k_7_r56 .* BMAL1nuc .* PERCRYnuc));
	k_7_r55 = (1.0 .* 1000.0 ./ 3600.0);
	J_r55 = ((k_8_r55 .* PERCRY_CLOCKBMAL1) - (k_7_r55 .* BMAL1nuc .* PERCRYnuc));
	J_r54 = ((nu_dPCN .* PERCRYnuc_p ./ (K_d_r54 + PERCRYnuc_p)) + (k_dn_r54 .* PERCRYnuc_p));
	J_r53 = (k_dn_r53 .* PERCRYnuc);
	J_r52 = ((V_3PC .* PERCRYnuc ./ (K_p_r52 + PERCRYnuc)) - (V_4PC .* PERCRYnuc_p ./ (K_dp_r52 + PERCRYnuc_p)));
	J_r51 = ( - (k_4_r51 .* PERCRY) + (k_3_r51 .* PER .* CRY));
	K_d_r40 = (0.3 ./ 1000.0);
	LumpedJ_r50 = (((k_1 .* PERCRY) - (k_2 .* PERCRYnuc)) .* Size_Cyto ./ KMOLE);
	J_r49 = ((nu_sC .* (BMAL1nuc ^ n_r49) ./ ((K_AC ^ n_r49) + (BMAL1nuc ^ n_r49))) + (nu_MRTF_C .* MRTFHillFcn) + nu_YTC.*YAPTAZHillFcn);
	J_r48 = ((nu_mC .* CRY_mRNA ./ (K_mC + CRY_mRNA)) + (k_dmc .* CRY_mRNA));
	J_r47 = ((nu_mP .* PER_mRNA ./ (K_mP + PER_mRNA)) + (k_dmp .* PER_mRNA));
	J_r46 = ((nu_sP .* (BMAL1nuc ^ n_r46) ./ ((K_AP ^ n_r46) + (BMAL1nuc ^ n_r46))) + (nu_MRTF_P .* MRTFHillFcn) + nu_YTP.*YAPTAZHillFcn);
	J_r45 = ((k_4_r45 .* PERCRY) - (k_3_r45 .* PER .* CRY));
	J_r43 = ((k_4_r43 .* PERCRY) - (k_3_r43 .* PER .* CRY));
	J_r41 = (k_dn_r41 .* PER);
	J_r40 = ((nu_dPC .* PER_p ./ (K_d_r40 + PER_p)) + (k_dn_r40 .* PER_p));
	% Rates
	dydt = [
		 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r59);    % rate for MRTF
		((KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r50 ./ Size_NM) - J_r52 - J_r53 + J_r55);    % rate for PERCRYnuc
		J_r15;    % rate for LaminA
		J_r9;    % rate for MyoA
		((KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r18 ./ Size_NM) - J_r20 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r21 ./ Size_NM) - J_r29);    % rate for BMAL1
		((KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r12) - J_r0 + J_r2);    % rate for Fakp
		0.0;    % rate for TEAD
		(J_r57 - J_r58);    % rate for PERCRY_CLOCKBMAL1
		 - J_r11;    % rate for Gactin
		( - (KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r12) + J_r0 - J_r2);    % rate for Fak
		( - (KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r4) + J_r6);    % rate for ROCK
		(J_r27 - J_r28);    % rate for BMAL1_mRNA
		J_r10;    % rate for CofilinNP
		((KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r17 ./ Size_NM) - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r22 ./ Size_NM) - J_r24);    % rate for REV
		( - J_r19 + (KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r21 ./ Size_NM) - J_r31 + J_r56);    % rate for BMAL1nuc
		J_r14;    % rate for NPCA
		(J_r46 - J_r47);    % rate for PER_mRNA
		0.0;    % rate for Positionboolean
		 - J_r8;    % rate for LIMK
		((KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r3) - (KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r1));    % rate for RhoAGDP
		(KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* J_r59);    % rate for MRTFnuc
		J_r11;    % rate for Fcyto
		(J_r19 - J_r32);    % rate for BMAL1nuc_p
		 - J_r9;    % rate for Myo
		(J_r33 - J_r36);    % rate for PERCRY_p
		 - J_r13;    % rate for YAPTAZP
		(J_r13 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r16));    % rate for YAPTAZN
		(J_r23 - J_r25);    % rate for REV_mRNA
		0.0;    % rate for SRF
		( - J_r34 - J_r39 + (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r44 ./ Size_NM) + J_r45);    % rate for CRY
		 - J_r15;    % rate for LaminAp
		(J_r52 - J_r54);    % rate for PERCRYnuc_p
		(KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* J_r16);    % rate for YAPTAZnuc
		(J_r35 - J_r40);    % rate for PER_p
		((KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r5) - J_r7);    % rate for mDiaA
		J_r8;    % rate for LIMKA
		( - (KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r5) + J_r7);    % rate for mDia
		( - J_r33 - J_r37 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r50 ./ Size_NM) + J_r51);    % rate for PERCRY
		( - J_r48 + J_r49);    % rate for CRY_mRNA
		((KFlux_PM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* J_r4) - J_r6);    % rate for ROCKA
		( - J_r3 + J_r1);    % rate for RhoAGTP_MEM
		( - J_r35 - J_r41 + (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r42 ./ Size_NM) + J_r43);    % rate for PER
		((KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r22 ./ Size_NM) - J_r26);    % rate for REVnuc
		 - J_r10;    % rate for CofilinP
		(J_r34 - J_r38);    % rate for CRY_p
		 - J_r14;    % rate for NPC
		(J_r20 - J_r30);    % rate for BMAL1_p
	];
end
