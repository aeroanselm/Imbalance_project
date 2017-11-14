clear all
close all
clc


% Defined Constants
g0 = 9.81;


% PAYLOAD FAIRING
m_fair = 2675;
m_pay = 10000;
m0 = 780000;


% CRYOGENIC MAIN CORE STAGE (EPC)
m_dry_EPC = 14700;
m_prop_EPC = 170000;                       % This value has to be corrected (estimated)
t_burn_EPC = 540;
Isp_EPC = 310;                             % This value has to be corrected (estimated)
Isp_EPC_vac = 432;
T_EPC = 960000;
T_EPC_vac = 1390000;
m_dot_EPC = T_EPC_vac/Isp_EPC_vac/g0;
Isp_EPC = T_EPC/m_dot_EPC/g0;              % Corrected value
m_prop_EPC = m_dot_EPC*t_burn_EPC;

%m_LOX = 150000;
%m_LH2 = 25000;


% SOLID ROCKET BOOSTER (EAP)
t_burn_EAP = 135;
m_prop_EAP = 237800;
Isp_EAP = 262;
T_EAP = 5400000;
m_dot_EAP_single = T_EAP/g0/Isp_EAP;
T_EAP_both = T_EAP*2;
m_dot_EAP = m_dot_EAP_single*2;


% CRYOGENIC UPPER STAGE (ESC-A)
m_dry_ESC = 4540;
t_burn_ESC = 945;
T_ESC = 67000;
Isp_ESC = 446;
m_dot_ESC = T_ESC/Isp_ESC/g0;
m_prop_ESC = m_dot_ESC*t_burn_ESC;


%% Vertical launch without drag. Thrust+gravity

T_par = T_EAP_both + T_EPC;

Isp_par = T_par/((m_dot_EAP + m_dot_EPC)*g0);
t_mission = 30;
mf_30s = m0 - m_dot_EAP*t_mission - m_dot_EPC*t_mission;
DV_paraller = -Isp_par*g0*log(mf_30s/m0) - g0*t_mission;

ve = T_par/(m_dot_EAP + m_dot_EPC);
MR = mf_30s/m0;
h = -ve*t_mission*(MR*log(MR)/(MR-1)-1) - 0.5*g0*t_mission;
