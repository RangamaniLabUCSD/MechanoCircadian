function [T,Y] = CircadianOnly_full(timeSpan, circad_param, mechanoVals)
% Implementation of Circadian model from Leloup and Goldbeter 2003
% Developed in VCell, output to this MATLAB function
%
% input:
%     timeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
%     circad_param is a normalized vector that multiplies each parameter in
%     this model
%     mechanoVals is a two element vector with the YAPTAZ nuclear
%     concentration (first el) and the MRTF nuclear conc (second el)
%
% output:
%     T is the vector of times
%     Y is the vector of state variables

%
% Default Initial Conditions
%
yinit = [
	0.002;		% yinit(1) is the initial condition for 'PERCRYnuc'
	0.005;		% yinit(2) is the initial condition for 'BMAL1'
	0.0;		% yinit(3) is the initial condition for 'PERCRY_CLOCKBMAL1'
	0.002;		% yinit(4) is the initial condition for 'BMAL1_mRNA'
	0.0;		% yinit(5) is the initial condition for 'REV'
	0.005;		% yinit(6) is the initial condition for 'BMAL1nuc'
	0.002;		% yinit(7) is the initial condition for 'PER_mRNA'
	0.0;		% yinit(8) is the initial condition for 'MRTFnuc'
	0.0;		% yinit(9) is the initial condition for 'BMAL1nuc_p'
	0.0;		% yinit(10) is the initial condition for 'PERCRY_p'
	0.002;		% yinit(11) is the initial condition for 'REV_mRNA'
	0.002;		% yinit(12) is the initial condition for 'CRY'
	0.0;		% yinit(13) is the initial condition for 'PERCRYnuc_p'
	0.7;		% yinit(14) is the initial condition for 'YAPTAZnuc'
	0.0;		% yinit(15) is the initial condition for 'PER_p'
	0.002;		% yinit(16) is the initial condition for 'PERCRY'
	0.002;		% yinit(17) is the initial condition for 'CRY_mRNA'
	0.002;		% yinit(18) is the initial condition for 'PER'
	0.0;		% yinit(19) is the initial condition for 'REVnuc'
	0.0;		% yinit(20) is the initial condition for 'CRY_p'
	0.0;		% yinit(21) is the initial condition for 'BMAL1_p'
];

yinit(14) = mechanoVals(1); % value of YAPTAZNuc
yinit(8) = mechanoVals(2); % value of MRTFNuc

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
	0.04;		% param(5) is 'kturnover'
	0.002;		% param(6) is 'PER_init_uM'
	1.0;		% param(7) is 'netValence_r50'
	0.0;		% param(8) is 'CRY_p_init_uM'
	0.8;		% param(9) is 'kdrock'
	0.002;		% param(10) is 'PERCRY_init_uM'
	0.0;		% param(11) is 'BMAL1_p_init_uM'
	7.6E-4;		% param(12) is 'kCY'
	1.0;		% param(13) is 'netValence_r44'
	1.0;		% param(14) is 'netValence_r42'
	0.56;		% param(15) is 'kCN'
	0.379;		% param(16) is 'ksf'
	4.0;		% param(17) is 'kfc1'
	0.067;		% param(18) is 'kdmy'
	0.648;		% param(19) is 'krp'
	0.0;		% param(20) is 'PER_p_init_uM'
	0.4;		% param(21) is 'kra'
	1.0;		% param(22) is 'netValence_r22'
	1.0;		% param(23) is 'netValence_r21'
	77.56;		% param(24) is 'gamma'
	1.0;		% param(25) is 'netValence_r18'
	1.0;		% param(26) is 'netValence_r17'
	0.002;		% param(27) is 'PERCRYnuc_init_uM'
	1.8;		% param(28) is 'SAV'
	1000.0;		% param(29) is 'K_millivolts_per_volt'
	1.0E-9;		% param(30) is 'mlabfix_K_GHK_'
	1.0;		% param(31) is 'q_r49'
	2.0;		% param(32) is 'q_r46'
	0.0;		% param(33) is 'MRTFnuc_init_uM'
	9.64853321E-5;		% param(34) is 'mlabfix_F_nmol_'
	393.0;		% param(35) is 'Size_NM'
	2.0;		% param(36) is 'q_r27'
	0.03;		% param(37) is 'kmr'
	0.002;		% param(38) is 'kmp'
	0.0;		% param(39) is 'REV_init_uM'
	120.0;		% param(40) is 'kly'
	0.07;		% param(41) is 'klr'
	16.0;		% param(42) is 'kll'
	300.0;		% param(43) is 'mlabfix_T_'
	2300.0;		% param(44) is 'Size_Cyto'
	550.0;		% param(45) is 'Size_Nuc'
	100.0;		% param(46) is 'Clamin'
	3.5;		% param(47) is 'kdep'
	2.0;		% param(48) is 'n_r49'
	2.0;		% param(49) is 'n_r46'
	0.002;		% param(50) is 'BMAL1_mRNA_init_uM'
	8314.46261815;		% param(51) is 'mlabfix_R_'
	0.015;		% param(52) is 'kf'
	0.0;		% param(53) is 'PERCRY_p_init_uM'
	0.7;		% param(54) is 'YAPTAZnuc_init_uM'
	0.002;		% param(55) is 'REV_mRNA_init_uM'
	55.49;		% param(56) is 'tau'
	0.0;		% param(57) is 'BMAL1nuc_p_init_uM'
	0.002;		% param(58) is 'CRY_mRNA_init_uM'
	0.005;		% param(59) is 'BMAL1_init_uM'
	8.7;		% param(60) is 'krNPC'
	0.0;		% param(61) is 'PERCRYnuc_p_init_uM'
	0.14;		% param(62) is 'kNC'
	0.165;		% param(63) is 'mDiaB'
	0.0;		% param(64) is 'REVnuc_init_uM'
	0.625;		% param(65) is 'kdp'
	0.0168;		% param(66) is 'kfkp'
	2.0;		% param(67) is 'kdl'
	0.035;		% param(68) is 'kdf'
	25.0;		% param(69) is 'Emol'
	2.8E-7;		% param(70) is 'kfNPC'
	3.141592653589793;		% param(71) is 'mlabfix_PI_'
	50.0;		% param(72) is 'alpha'
	0.3;		% param(73) is 'ROCKB'
	0.002;		% param(74) is 'PER_mRNA_init_uM'
	602.2;		% param(75) is 'unitconversionfactor'
	2.0;		% param(76) is 'm'
	0.005;		% param(77) is 'kdmdia'
	2.0;		% param(78) is 'h'
	0.001660538783162726;		% param(79) is 'KMOLE'
	6.02214179E11;		% param(80) is 'mlabfix_N_pmol_'
	16.0;		% param(81) is 'k11'
	3.25;		% param(82) is 'C'
	2.0;		% param(83) is 'K_YTB'
	36.0;		% param(84) is 'epsilon'
	0.005;		% param(85) is 'BMAL1nuc_init_uM'
	0.0;		% param(86) is 'PERCRY_CLOCKBMAL1_init_uM'
	1.0;		% param(87) is 'K_MRTF_r49'
];

%
% invoke the integrator
%
[T,Y] = ode15s(@f,timeSpan,yinit,odeset(),param,yinit,circad_param);

end


% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0, circad_param)
	% State Variables
	PERCRYnuc = y(1);
	BMAL1 = y(2);
	PERCRY_CLOCKBMAL1 = y(3);
	BMAL1_mRNA = y(4);
	REV = y(5);
	BMAL1nuc = y(6);
	PER_mRNA = y(7);
	MRTFnuc = y(8);
	BMAL1nuc_p = y(9);
	PERCRY_p = y(10);
	REV_mRNA = y(11);
	CRY = y(12);
	PERCRYnuc_p = y(13);
	YAPTAZnuc = y(14);
	PER_p = y(15);
	PERCRY = y(16);
	CRY_mRNA = y(17);
	PER = y(18);
	REVnuc = y(19);
	CRY_p = y(20);
	BMAL1_p = y(21);

	% Constants
	q_r49 = p(31);
	Size_NM = p(35);
	q_r27 = p(36);
	Size_Cyto = p(44);
	Size_Nuc = p(45);
	n_r49 = p(48);
	n_r46 = p(49);
	m = p(76);
	h = p(78);
	KMOLE = p(79);
	K_YTB = p(83);
	K_MRTF_r49 = p(87);
    
    couplingRef = 0.5 / (1000.0 * 3600.0);
    nu_YTB = couplingRef *    circad_param(1);
    nu_YTP = couplingRef *    circad_param(2);
    nu_YTC = couplingRef *    circad_param(3);
    nu_YTR = couplingRef *    circad_param(4);
    nu_MRTF_B = couplingRef * circad_param(5);
	nu_MRTF_P = couplingRef * circad_param(6);
    nu_MRTF_C = couplingRef * circad_param(7);
    nu_MRTF_R = couplingRef * circad_param(8);
    YAPTAZHillFcn = ((YAPTAZnuc ^ q_r27)) ./ ((K_YTB ^ q_r27) + (YAPTAZnuc ^ q_r27));
    MRTFHillFcn = (MRTFnuc ^ q_r49) ./ ((K_MRTF_r49 ^ q_r49) + (MRTFnuc ^ q_r49));

	nu_dRC = (4.4 .* 0.001 ./ 3600.0)        * circad_param(9);
	K_d_r36 = (0.3 ./ 1000.0)                * circad_param(10); 
	K_d_r32 = (0.3 .* 0.001)                 * circad_param(11);
	k_sC = (3.2 ./ 3600.0)                   * circad_param(12);
	K_d_r30 = (0.3 .* 0.001)                 * circad_param(13);
	k_sP = (1.2 ./ 3600.0)                   * circad_param(14);
	k_9 = (0.8 ./ 3600.0)                    * circad_param(15);
	k_6 = (0.8 ./ 3600.0)                    * circad_param(16);
	k_5 = (0.8 ./ 3600.0)                    * circad_param(17);
	k_dnC = (0.02 ./ 3600.0)                 * circad_param(18);
	K_d_r38 = (0.3 ./ 1000.0)                * circad_param(19);
	k_dn_r38 = (0.02 ./ 3600.0)              * circad_param(20);  
	nu_dCC = (1.4 ./ (1000.0 .* 3600.0))     * circad_param(21);
	k_dn_r37 = (0.02 ./ 3600.0)              * circad_param(22);
	k_2 = (0.4 ./ 3600.0)                    * circad_param(23); 
	nu_dPCC = (1.4 ./ (1000.0 .* 3600.0))    * circad_param(24); 
	k_dn_r36 = (0.02 ./ 3600.0)              * circad_param(25);
	k_1 = (0.8 ./ 3600.0)                    * circad_param(26);
	V_1P = (9.6 ./ (1000.0 .* 3600.0))       * circad_param(27);
	K_dp_r35 = (0.1 ./ 1000.0)               * circad_param(28);
	K_p_r35 = (1.006 ./ 1000.0)              * circad_param(29);
	V_2P = (0.6 ./ (1000.0 .* 3600.0))       * circad_param(30);
	V_2C = (0.2 ./ (1000.0 .* 3600.0))       * circad_param(31);
	V_1C = (1.2 ./ (1000.0 .* 3600.0))       * circad_param(33);
	K_p_r34 = (1.006 ./ 1000.0)              * circad_param(34);
	K_dp_r34 = (0.1 ./ 1000.0)               * circad_param(35);
	K_dp_r33 = (0.1 ./ 1000.0)               * circad_param(36);
	K_p_r33 = (1.006 ./ 1000.0)              * circad_param(37);
	V_2PC = (0.2 ./ (1000.0 .* 3600.0))      * circad_param(38);
	V_1PC = (2.4 ./ (1000.0 .* 3600.0))      * circad_param(39);
	nu_dBN = (3.0 .* 0.001 ./ 3600.0)        * circad_param(40);
	k_dn_r32 = (0.02 ./ 3600.0)              * circad_param(41);
	k_dn_r31 = (0.02 ./ 3600.0)              * circad_param(42);
	k_dn_r30 = (0.02 ./ 3600.0)              * circad_param(43);
	nu_dBC = (3.0 .* 0.001 ./ 3600.0)        * circad_param(44);
	K_d_r26 = (0.3 .* 0.001)                 * circad_param(45);
	K_d_r24 = (0.3 .* 0.001)                 * circad_param(46);
	k_dn_r29 = (0.02 ./ 3600.0)              * circad_param(47);
	k_dmB = (0.02 ./ 3600.0)                 * circad_param(48);
	nu_mB = (1.3 .* 0.001 ./ 3600.0)         * circad_param(49);
	K_mB = (0.4 .* 0.001)                    * circad_param(50);
	nu_sB = (1.8 .* 0.001 ./ 3600.0)         * circad_param(51);
	K_IB = (2.2 .* 0.001)                    * circad_param(52);
	k_dn_r26 = (0.02 ./ 3600.0)              * circad_param(53);
	nu_dRN = (0.8 .* 0.001 ./ 3600.0)        * circad_param(54);
	k_dmr = (0.02 ./ 3600.0)                 * circad_param(55);
	nu_mR = (1.6 .* 0.001 ./ 3600.0)         * circad_param(56);
	K_mR = (0.4 .* 0.001)                    * circad_param(57);
	k_dn_r24 = (0.02 ./ 3600.0)              * circad_param(58);
	nu_sR = (1.6 .* 0.001 ./ 3600.0)         * circad_param(59);
	K_AR = (0.6 .* 0.001)                    * circad_param(60);
	V_2B = (0.2 .* 0.001 ./ 3600.0)          * circad_param(61);
	V_1B = (1.4 .* 0.001 ./ 3600.0)          * circad_param(62);
	K_dp_r20 = (0.1 .* 0.001)                * circad_param(63);
	K_p_r20 = (1.006 .* 0.001)               * circad_param(64);
	nu_dPC = (3.4 ./ (1000.0 .* 3600.0))     * circad_param(65);
	k_10 = (0.4 ./ 3600.0)                   * circad_param(66);
	K_p_r52 = (1.006 ./ 1000.0)              * circad_param(67);
	V_3B = (1.4 .* 0.001 ./ 3600.0)          * circad_param(68);
	K_p_r19 = (1.006 .* 0.001)               * circad_param(69);
	K_dp_r19 = (0.1 .* 0.001)                * circad_param(70);
	V_4B = (0.4 .* 0.001 ./ 3600.0)          * circad_param(71);
	k_sB = (0.32 ./ 3600.0)                  * circad_param(72);
	k_sR = (1.7 ./ 3600.0)                   * circad_param(73);
	K_AP = (0.6 ./ 1000.0)                   * circad_param(74);
	K_AC = (0.6 ./ 1000.0)                   * circad_param(75);
	k_4_r51 = (0.4 ./ 3600.0)                * circad_param(76);
	nu_sP = (2.4 ./ (1000.0 .* 3600.0))      * circad_param(77);
	k_4_r45 = (0.4 ./ 3600.0)                * circad_param(78);
	k_4_r43 = (0.4 ./ 3600.0)                * circad_param(79);
	nu_sC = (2.2 ./ (1000.0 .* 3600.0))      * circad_param(80);
	k_3_r51 = (0.8 .* 1000.0 ./ 3600.0)      * circad_param(81);
	nu_dPCN = (1.4 ./ (1000.0 .* 3600.0))    * circad_param(82);
	k_3_r45 = (0.8 .* 1000.0 ./ 3600.0)      * circad_param(83);
	k_3_r43 = (0.8 .* 1000.0 ./ 3600.0)      * circad_param(84);
	nu_dIN = (1.6 ./ (1000.0 .* 3600.0))     * circad_param(85);
	K_mP = (0.3 ./ 1000.0)                   * circad_param(86);
	K_mC = (0.4 ./ 1000.0)                   * circad_param(87);
	nu_mP = (2.2 ./ (1000.0 .* 3600.0))      * circad_param(88);
	nu_mC = (2.0 ./ (1000.0 .* 3600.0))      * circad_param(89);
	k_dn_r58 = (0.02 ./ 3600.0)              * circad_param(90);
	k_dn_r54 = (0.02 ./ 3600.0)              * circad_param(91);
	k_dn_r53 = (0.02 ./ 3600.0)              * circad_param(92);
	k_dn_r41 = (0.02 ./ 3600.0)              * circad_param(93);
	k_dn_r40 = (0.02 ./ 3600.0)              * circad_param(94);
	K_dp_r52 = (0.1 ./ 1000.0)               * circad_param(95);
	V_4PC = (0.2 ./ (1000.0 .* 3600.0))      * circad_param(96);
	k_dmp = (0.02 ./ 3600.0)                 * circad_param(97);
	k_dmc = (0.02 ./ 3600.0)                 * circad_param(98);
	k_8_r57 = (0.2 ./ 3600.0)                * circad_param(99);
	K_d_r58 = (0.3 ./ 1000.0)                * circad_param(100);
	k_8_r56 = (0.2 ./ 3600.0)                * circad_param(101);
	k_8_r55 = (0.2 ./ 3600.0)                * circad_param(102);
	K_d_r54 = (0.3 ./ 1000.0)                * circad_param(103);
	V_3PC = (2.4 ./ (1000.0 .* 3600.0))      * circad_param(104);
	k_7_r57 = (1.0 .* 1000.0 ./ 3600.0)      * circad_param(105);
	k_7_r56 = (1.0 .* 1000.0 ./ 3600.0)      * circad_param(106);
	k_7_r55 = (1.0 .* 1000.0 ./ 3600.0)      * circad_param(107);
	K_d_r40 = (0.3 ./ 1000.0)                * circad_param(108);


	% Functions
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 .* (KMOLE ^ 1.0));
	LumpedJ_r44 = (k_sC .* CRY_mRNA .* Size_Cyto ./ KMOLE);
	LumpedJ_r42 = ((k_sP .* PER_mRNA) .* Size_Cyto ./ KMOLE);
	J_r39 = (k_dnC .* CRY);
	J_r38 = ((nu_dCC .* CRY_p ./ (K_d_r38 + CRY_p)) + (k_dn_r38 .* CRY_p));
	J_r37 = (k_dn_r37 .* PERCRY);
	J_r36 = ((nu_dPCC .* PERCRY_p ./ (K_d_r36 + PERCRY_p)) + (k_dn_r36 .* PERCRY_p));
	J_r35 = ((V_1P .* PER ./ (K_p_r35 + PER)) - (V_2P .* PER_p ./ (K_dp_r35 + PER_p)));
	J_r34 = ((V_1C .* CRY ./ (K_p_r34 + CRY)) - (V_2C .* CRY_p ./ (K_dp_r34 + CRY_p)));
	J_r33 = ((V_1PC .* PERCRY ./ (K_p_r33 + PERCRY)) - (V_2PC .* PERCRY_p ./ (K_dp_r33 + PERCRY_p)));
	J_r32 = ((nu_dBN .* BMAL1nuc_p ./ (K_d_r32 + BMAL1nuc_p)) + (k_dn_r32 .* BMAL1nuc_p));
	J_r31 = (k_dn_r31 .* BMAL1nuc);
	J_r30 = ((nu_dBC .* BMAL1_p ./ (K_d_r30 + BMAL1_p)) + (k_dn_r30 .* BMAL1_p));
	J_r29 = (k_dn_r29 .* BMAL1);
	J_r28 = (((nu_mB .* BMAL1_mRNA) ./ (K_mB + BMAL1_mRNA)) + (k_dmB .* BMAL1_mRNA));
	J_r27 = (((nu_sB .* (K_IB ^ m)) ./ ((K_IB ^ m) + (REVnuc ^ m))) + (nu_YTB .* YAPTAZHillFcn) + nu_MRTF_B.*MRTFHillFcn);
	J_r26 = (((nu_dRN .* REVnuc) ./ (K_d_r26 + REVnuc)) - (k_dn_r26 .* REVnuc));
	J_r25 = (((nu_mR .* REV_mRNA) ./ (K_mR + REV_mRNA)) + (k_dmr .* REV_mRNA)) + nu_YTR.*YAPTAZHillFcn + nu_MRTF_R.*MRTFHillFcn;
	J_r24 = (((nu_dRC .* REV) ./ (K_d_r24 + REV)) + (k_dn_r24 .* REV));
	J_r23 = ((nu_sR .* (BMAL1nuc ^ h)) ./ ((K_AR ^ h) + (BMAL1nuc ^ h)));
	J_r20 = ((V_1B .* BMAL1 ./ (K_p_r20 + BMAL1)) - (V_2B .* BMAL1_p ./ (K_dp_r20 + BMAL1_p)));
	LumpedJ_r22 = (((k_9 .* REV) - (k_10 .* REVnuc)) .* Size_Cyto ./ KMOLE);
	LumpedJ_r21 = (((k_5 .* BMAL1) - (k_6 .* BMAL1nuc)) .* Size_Cyto ./ KMOLE);
	J_r19 = ((V_3B .* BMAL1nuc ./ (K_p_r19 + BMAL1nuc)) - (V_4B .* BMAL1nuc_p ./ (K_dp_r19 + BMAL1nuc_p)));
	KFlux_NM_Nuc = (Size_NM ./ Size_Nuc);
	LumpedJ_r18 = (k_sB .* BMAL1_mRNA .* Size_Cyto ./ KMOLE);
	LumpedJ_r17 = (k_sR .* REV_mRNA .* Size_Cyto ./ KMOLE);
	KFlux_NM_Cyto = (Size_NM ./ Size_Cyto);
	J_r58 = ((nu_dIN .* PERCRY_CLOCKBMAL1 ./ (K_d_r58 + PERCRY_CLOCKBMAL1)) + (k_dn_r58 .* PERCRY_CLOCKBMAL1));
	J_r57 = ( - (k_8_r57 .* PERCRY_CLOCKBMAL1) + (k_7_r57 .* BMAL1nuc .* PERCRYnuc));
	J_r56 = ((k_8_r56 .* PERCRY_CLOCKBMAL1) - (k_7_r56 .* BMAL1nuc .* PERCRYnuc));
	J_r55 = ((k_8_r55 .* PERCRY_CLOCKBMAL1) - (k_7_r55 .* BMAL1nuc .* PERCRYnuc));
	J_r54 = ((nu_dPCN .* PERCRYnuc_p ./ (K_d_r54 + PERCRYnuc_p)) + (k_dn_r54 .* PERCRYnuc_p));
	J_r53 = (k_dn_r53 .* PERCRYnuc);
	J_r52 = ((V_3PC .* PERCRYnuc ./ (K_p_r52 + PERCRYnuc)) - (V_4PC .* PERCRYnuc_p ./ (K_dp_r52 + PERCRYnuc_p)));
	J_r51 = ( - (k_4_r51 .* PERCRY) + (k_3_r51 .* PER .* CRY));
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
		((KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r50 ./ Size_NM) - J_r52 - J_r53 + J_r55);    % rate for PERCRYnuc
		((KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r18 ./ Size_NM) - J_r20 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r21 ./ Size_NM) - J_r29);    % rate for BMAL1
		(J_r57 - J_r58);    % rate for PERCRY_CLOCKBMAL1
		(J_r27 - J_r28);    % rate for BMAL1_mRNA
		((KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r17 ./ Size_NM) - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r22 ./ Size_NM) - J_r24);    % rate for REV
		( - J_r19 + (KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r21 ./ Size_NM) - J_r31 + J_r56);    % rate for BMAL1nuc
		(J_r46 - J_r47);    % rate for PER_mRNA
		0.0;    % rate for MRTFnuc
		(J_r19 - J_r32);    % rate for BMAL1nuc_p
		(J_r33 - J_r36);    % rate for PERCRY_p
		(J_r23 - J_r25);    % rate for REV_mRNA
		( - J_r34 - J_r39 + (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r44 ./ Size_NM) + J_r45);    % rate for CRY
		(J_r52 - J_r54);    % rate for PERCRYnuc_p
		0.0;    % rate for YAPTAZnuc
		(J_r35 - J_r40);    % rate for PER_p
		( - J_r33 - J_r37 - (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r50 ./ Size_NM) + J_r51);    % rate for PERCRY
		( - J_r48 + J_r49);    % rate for CRY_mRNA
		( - J_r35 - J_r41 + (KFlux_NM_Cyto .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r42 ./ Size_NM) + J_r43);    % rate for PER
		((KFlux_NM_Nuc .* UnitFactor_uM_um3_molecules_neg_1 .* LumpedJ_r22 ./ Size_NM) - J_r26);    % rate for REVnuc
		(J_r34 - J_r38);    % rate for CRY_p
		(J_r20 - J_r30);    % rate for BMAL1_p
	];
end
