function [SSVar, tauVals] = MechanoSS(stiffnessParam, inhibVec, pSol, varargin)
    % function returning the steady state state variable values for YAP/TAZ
    % and MRTF mechanotransduction model (developed from VCell)

    % inputs:
    %     stiffnessParam is a 2 element vector with the stiffness in kPa (first
    %       element) and the time scale of stiffness increase in s (second el)
    %       only stiffnessParam(1) is used here (stiffness is assumed to be
    %       equal to this value at SS)
    %
    %     inhibVec: vector of inhibition parameters (length=9)
    %       1: actin polym inhibition: factor multiplying kra
    %       2: ROCK inhibition: factor multiplying param 55, 69 (epsilon and tau - ROCK mediated catalysis)
    %       3: MRTF-Circadian coupling inhibition (not used here)
    %       4: YAP overexpression - fold expression of 5SA-YAP (compared to normal YAP expression)
    %       5: CytD concentration (in micromolar)
    %       6: LATS factor - factor multiplying kNC (rate of YAP phosphorylation)
    %       7: blebbistatin - factor multiplying param(46) (rate of stress fiber dependent YAPTAZ
    %           dephos) and param(86) (rate of stress fiber-dependent nuclear pore opening)
    %       8: cell contact area (in microns squared, control area is 3000)
    %       9: lamin A mutation - factor multiplying lamin A phos rate (krl)
    %
    %     pSol: parameter solution vector (see full list in
    %     MechanoCircadianModel function)
    %
    %     popVar (optional) - either length of 105 (total # of parameters
    %     to vary) or scalar corresponding to population variance
    %
    %
    % Outputs:
    %   SSVar: SS values for all state variables
    %   tauVals: characteristic time scale computed for each variable

    if isempty(varargin)
        popVar = 0;
    elseif length(varargin)==1
        popVar = varargin{1};
    end

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
        0.001;		% param(52) is 'krl'
        0.0;		% param(53) is 'Kr_r12'
        100.0;		% param(54) is 'Clamin'
        36.0;		% param(55) is 'epsilon'
        0.2;		% param(56) is 'YAPTAZP_init_uM'
        17.9;		% param(57) is 'Fcyto_init_uM'
        0.3;		% param(58) is 'ROCKB'
        482.4;		% param(59) is 'Gactin_init_uM'
        1080.0;		% param(60) is 'Size_PM'
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
        0.0;		% param(80) is 'Kr_r3
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
    param(104) = stiffnessParam(2); % unused
    param(78) = pSol(19); % substrate stiffness sensitivity (C)
    
    % load inhibition parameters 
    % [actin polym (kra), ROCK, MRTF, YAP overexpression, cytD conc, LATS factor,...
    % blebbistatin conc, contact area, Lamin A phos]
    param(3) = param(3)*inhibVec(1); % inhibition of actin polymerization
    param([55,69]) = param([55,69])*inhibVec(2); % inhibit ROCK-mediated catalysis
    param(112) = inhibVec(4); % fold 5SA-YAP overexpression
    param(105) = inhibVec(5); % cytD conc
    param(2) = inhibVec(6)*param(2); % LATS factor (for changes in cell density)
    param(86) = inhibVec(7)*param(86); % inhibit actin-myosin interactions (blebbistatin) - YAPTAZ dephos
    param(46) = inhibVec(7)*param(46); % inhibit actin-myosin interactions (blebbistatin) - NPC opening
    param([5,97]) = inhibVec(8)*param([5,97])/3000; % account for contact area
    param(52) = inhibVec(9)*param(52); % inhibit phosph. of LaminA
    
    param(113) = pSol(30); %MRTFReleaseConst
    param(106) = 1e6; %MRTFTot
    param(107) = pSol(31); %Kinsolo_MRTF
    param(108) = pSol(32); %Kin2_MRTF
    param(109) = pSol(34); %Kcap (cytoD-dependent F actin capping rate)
    param(110) = pSol(33); %Kdim (cytoD-dependent G actin dimerization rate)

    % introduce variability in parameters (all rate constants and
    % concentrations except calibrated MRTF parameters, those are sampled
    % separately from their own distributions)
    if length(popVar) == 105
        varVecCur = popVar(1:105);
        varVecCur = varVecCur(:);
    else
        varVecCur = sqrt(popVar).*randn(size(param([1:104,106])));
    end
    varVecCur([33:36,47:49,51,64,66,68,72,85,90:92,94,99:100,102:103]) = 0; % exclude constants and exponents
    param([1:104,106]) = param([1:104,106]) .* exp(varVecCur);
    p = param;

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
	% Size_ECM = p(26);
	% Voltage_PM = p(27);
	kmr = p(28);
	gamma = p(29);
	kmp = p(30);
	CofilinNP_init_uM = p(31);
	% k11 = p(32);
	% netValence_r5 = p(33);
	% netValence_r4 = p(34);
	% netValence_r3 = p(35);
	% netValence_r1 = p(36);
	kflaminA = p(37);
	% kly = p(38);
	Fakp_init_uM = p(39);
	klr = p(40);
	% kll = p(41);
	LaminA_init_molecules_um_2 = p(42);
	mDiaB = p(43);
	% Voltage_NM = p(44);
	RhoAGTP_MEM_init_molecules_um_2 = p(45);
	kfNPC = p(46);
	% mlabfix_F_nmol_ = p(47);
	unitconversionfactor = p(48);
	% mlabfix_T_ = p(49);
	mDia_init_uM = p(50);
	% K_millivolts_per_volt = p(51);
	krl = p(52);
	% Kr_r12 = p(53);
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
	% mlabfix_PI_ = p(64);
	propStiff = p(65);
	% mlabfix_F_ = p(66);
	kinSolo2 = p(67);
	% mlabfix_R_ = p(68);
	tau = p(69);
	kout2 = p(70);
	ROCK_init_uM = p(71);
	% mlabfix_K_GHK_ = p(72);
	kdep = p(73);
	% Size_NM = p(74);
	% Kr_r7 = p(75);
	% Kr_r6 = p(76);
	% Kr_r5 = p(77);
	C = p(78);
	% Kr_r4 = p(79);
	% Kr_r3 = p(80);
	% Kr_r2 = p(81);
	% Kr_r1 = p(82);
	% Kr_r0 = p(83);
	NPC_init_molecules_um_2 = p(84);
	% mlabfix_N_pmol_ = p(85);
	kCY = p(86);
	CofilinP_init_uM = p(87);
	LaminAp_init_molecules_um_2 = p(88);
	kCN = p(89);
	% netValence_r16 = p(90);
	% netValence_r15 = p(91);
	% netValence_r14 = p(92);
	LIMK_init_uM = p(93);
	% netValence_r12 = p(94);
	kturnover = p(95);
	Fak_init_uM = p(96);
	ksf = p(97);
	YAPTAZN_init_uM = p(98);
	n2 = p(99);
	n1 = p(100);
	NPCA_init_molecules_um_2 = p(101);
	% Positionboolean_init_molecules_um_2 = p(102);
	% KMOLE = p(103);
	KFlux_PM_Cyto = (Size_PM ./ Size_Cyto);
	% KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	UnitFactor_uM_um3_molecules_neg_1 = (1000000.0 ./ 6.02214179E8);

    cytoDConc = p(105);
    MRTFTot = p(106);
    Kinsolo_MRTF = p(107);
    Kin2_MRTF  = p(108);
    Kcap = p(109);
    Kdim = p(110);
    foldYAPOverexpress = p(112);
    MRTFReleaseConst = p(113);
	
    %SS calcs
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
    % GactinAlt = Fcyto*(kdep + kfc1*CofilinNP)/(kra*(alpha*smoothmDiaA+1));
    ActinTau = (kdep_tot + kra*(alpha*smoothmDiaA+1)*(1+cytoDConc/Kcap))^(-1);

    MyoTot = Myo_init_uM + MyoA_init_uM;
    MyoA = MyoTot*kmr*(epsilon*smoothROCKA+1)/(kdmy + kmr*(epsilon*smoothROCKA+1));
    Myo = MyoTot - MyoA;
    MyoTau = (kdmy + kmr*(epsilon*smoothROCKA+1))^(-1);

    Ecytosol = propStiff*(Fcyto^n1);
    LaminATot = LaminAp_init_molecules_um_2 + LaminA_init_molecules_um_2;
    LaminA = LaminATot*kflaminA*Ecytosol/(kflaminA*Ecytosol + krl*(Clamin+Ecytosol));
    LaminAp = LaminATot - LaminA;
    LaminATau = (kflaminA*Ecytosol/(Clamin+Ecytosol) + krl)^(-1);

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
    YAPTAZCytoTot = YAPTAZTot / (1 + (1-YAPTAZPFraction)*YAPTAZnucFraction/(1-YAPTAZnucFraction));
    YAPTAZP = YAPTAZCytoTot * YAPTAZPFraction;
    YAPTAZN = YAPTAZCytoTot - YAPTAZP;
    YAPTAZnuc = YAPTAZTot - YAPTAZCytoTot;
    if foldYAPOverexpress > 0 % overexpression of YAP mutant
        YAPTAZTot_5SA = YAPTAZTot*foldYAPOverexpress;
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

    %MRTF
    MRTFFactor = 1/(1+(Gactin/MRTFReleaseConst)^2);
    nucPerm_MRTF = Kinsolo_MRTF + Kin2_MRTF*NPCA;
    MRTFnuc = (MRTFTot*nucPerm_MRTF*MRTFFactor/CytoConvert) / (nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + 1);
    MRTFcyto = (MRTFTot - MRTFnuc*NucConvert)/CytoConvert;
    MRTFnucTau = (nucPerm_MRTF*MRTFFactor*NucConvert/CytoConvert + 1)^(-1);

	SSVar = [CofilinP; Fak; mDia; LaminA; Fcyto; RhoAGTP_MEM; mDiaA; NPCA;...
        Gactin; NPC; ROCKA; Myo; CofilinNP; LaminAp; YAPTAZnuc; Fakp;...
        YAPTAZP; YAPTAZN; RhoAGDP; LIMK; MyoA; ROCK; 1; LIMKA; MRTFnuc; MRTFcyto]';
    tauVals = [CofilinTau; FakTau; mDiaTau; LaminATau; ActinTau; RhoATau; mDiaTau; NPCTau;...
        ActinTau; NPCTau; ROCKTau; MyoTau; CofilinTau; LaminATau; YAPTAZnucTau; FakTau;...
        YAPTAZPTau; YAPTAZPTau; RhoATau; LIMKTau; MyoTau; ROCKTau; 1; LIMKTau; MRTFnucTau; MRTFnucTau]';
    
end