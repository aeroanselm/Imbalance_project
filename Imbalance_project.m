% Imbalance_project.m - Ariane 5 ECA ---Imbalance project---
%  
% TYPE: 
%   Main Script
%
% Notes for upgrading this script: 
%   Please do not upgrade directly this script. Make a your local copy and do
%   your modifications and upgrades and send me the upgraded script.
%   Contact the author for modifications.
%
% DESCRIPTION: 
%   This routine is used to simulate and analyze scenarios in
%   case of imbalance between the twin externaln rockets operating as boost 
%   stage relative to the Ariane 5 ECA launcher by ESA.         
%
% REFERENCES:
%   - Ariane 5 Users Manual October 2016
%   - http://www.b14643.de/Spacerockets_1/West_Europe/Ariane-5/Description/Frame.htm
%   - http://spaceflight101.com/ariane-5-va226/ariane-5-va226-launch-profile/
%   - Thrust profile of the EAP solid rocket boosters is taken from document
%     A5prop.pdf provided along with this script.
%
% FUTURE DEVELOPMENT:
%   Work in progress.
%
% ORIGINAL VERSION:
%   23/11/2017, ALESSANDRO MARIA MASSERINI
%				ALESSIO NEGRI
% --------------------------------------------------------------------------
% ***** BEGIN GPL LICENSE BLOCK *****
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ***** END GPL LICENCE BLOCK *****
% --------------------------------------------------------------------------
%
% AUTHORS:
%
%   Name: ALESSANDRO MARIA 
%   Surname: MASSERINI
%   ID number: 883251
%   Contact: alessandro.masserini@mail.polimi.it
%
%   Name: ALESSIO 
%   Surname: NEGRI
%   ID number: 874935
%   Contact: alessio1.negri@mail.polimi.it
%
%   Course: Space Engineering
%   Department: DAER
%   University: Politecnico di Milano
%   Class: Launch Systems 
%   Creation: 23/11/2017
%
% CHANGELOG: 
%   ALESSANDRO MARIA MASSERINI creation of this script - 23/11/2017
%
% --------------------------------------------------------------------------

clear 
close all
clc
%% Imbalance project
data
EAP
EPC
%% EPC - Mass flow rates of LOX and LH2 
of = 6;
mdot_LH2=EPC_mdot_p/(of+1);
mdot_LOX=of*mdot_LH2;

%% EAP - Propellant mass consumption at time instant (wrt EAP) 18 sec 71 sec
%  EAP start burning at T+7 seconds so the real instants of measurments are 11 and 64 seconds.  

T_boost = @(k) interp1(x,y,k);

EAP_mprop_discharged11 = EAP_mprop_int (EAP_mdot(0:18));
EAP_mprop_discharged64 = EAP_mprop_int (EAP_mdot(0:71));

%% EPC - Propellant mass consumption at time instant (wrt EAP) 18 sec 71 sec
%  EAP start burning at T+7 seconds so the real instants of measurments are 11 and 64 seconds.  

EPC_mprop_LH2_discharged11 = mdot_LH2*11;
EPC_mprop_LOX_discharged11 = mdot_LOX*11;
EPC_mprop_LH2_discharged64 = mdot_LH2*64;
EPC_mprop_LOX_discharged64 = mdot_LOX*64;


%% Aerodynamics 
% CALCULATIONS of areas

a_2 = d_EPC+2*d_EAP;                                                        % Major asis [m]
a = a_2/2;                                                                  % Semi-major axis [m]
b = d_ESCA/2;                                                               % Semi-minor axis [m]
A_nose = pi/4*d_c^2;                                                        % Area of the fictitious hemispherical nose [m^2]
A_EAP = pi/4*d_EAP^2;                                                       % EAP area [m^2]
A_EPC = pi/4*d_EPC^2;                                                       % EPC area [m^2]
A_ref = A_EPC + 2*A_EAP;                                                    % Reference area [m^2]
d_e_Vulcain = 2.10;                                                         % Vulcain exit nozzle diameter [m]
A_e_EPC = pi/4*d_e_Vulcain^2;                                               % Vulcain exit area [m^2]
d_e_EAP = 3.10;                                                             % EAP exit nozzle diameter [m]
A_e_EAP = pi/4*d_e_EAP^2;                                                   % EAP exit area [m^2]
A_e_bar = A_e_EPC +2*A_e_EAP;                                               % Sum of all the exit nozzle areas [m^2]
A_b = pi*a*b;                                                               % Elliptical base area [m^2]
A_f = pi/4*d_f^2;                                                           % Flare base area [m^2]

% C_D CALCULATIONS

load('c_D_vs_Ma')
c_D_data = c_D_vs_Ma(:,2);                                                  % c_D real data
Ma_data = c_D_vs_Ma(:,1);                                                   % Ma real data
Ma = linspace(0,6,1000);
c_D = spline(Ma_data,c_D_data,Ma);

figure()
plot(Ma,cd_Ma(Ma),'b')
grid on
ylim([0,1])
title('Ariane 5 ECA zero-lift drag coefficient')
xlabel('M')
ylabel('c_D_0')

% C_N & C_N_ALPHA CALCULATIONS

alpha_v = linspace(0,pi/2,1000);
c_N = @(alpha,phi) (a/b*cos(phi)^2 +...
    b/a*sin(phi)^2)*(abs(sin(2*alpha).*cos(alpha/2))+ ...
    1.3*l_B/(2*sqrt(a*b))*sin(alpha).^2);
c_N_alpha = @(alpha,phi) (a/b*cos(phi)^2 +...
    b/a*sin(phi)^2)*(abs(2*cos(2*alpha).*cos(alpha/2)-...
    0.5*sin(2*alpha).*sin(alpha./2)) +...
    1.3*l_B/(2*sqrt(a*b))*2*sin(alpha).*cos(alpha));

figure()
plot(alpha_v*180/pi,c_N(alpha_v,0),'b',...
     alpha_v*180/pi,c_N(alpha_v,5*pi/180),'r',...
     alpha_v*180/pi,c_N(alpha_v,10*pi/180),'g')
grid on
title('Normal force coefficient')
xlabel('\alpha [°]')
ylabel('c_N')
legend('\phi = 0°','\phi = 5°','\phi = 10°','location','northwest')

figure()
plot(alpha_v*180/pi,c_N_alpha(alpha_v,0),'b',...
     alpha_v*180/pi,c_N_alpha(alpha_v,5*pi/180),'r',...
     alpha_v*180/pi,c_N_alpha(alpha_v,10*pi/180),'g')
grid on
title('Normal force coefficient derived by \alpha')
xlabel('\alpha [°]')
ylabel('c_N_/_\alpha')
legend('\phi = 0°','\phi = 5°','\phi = 10°','location','northwest')

% X_CP CALCULATIONS

x_CP_mb = 2/3*l_N;
x_CP_EAP = 2/3*l_N_EAP*A_EAP/A_f + (1-A_EAP/A_f)*l_EAP - h_f*((d_m/d_EAP)^2 - 1)*A_EAP/A_f;
x_CP = (x_CP_mb*(l_B-l_EAP) + (x_CP_EAP + l_B-l_EAP)*2*l_EAP)/(l_B-l_EAP+2*l_EAP);

%% TVC

% Data at 18 s and 394.7 m, and 71 s and 14570 m after launch

CG = [16.23411 17.38154];                                                   % CG positions from the base [m]
CP_CG = l_B - x_CP - CG;                                                    % CP-CG distance [m]
T_EAP = [6485000 5242000];                                                  % EAP thrust [N]
T_EPC = [934500 1251000];                                                   % EPC thrust
q = [4686 34000];                                                           % Dynamic pressure [Pa]
y_bar = 2.7 + 1.525;                                                        % Point of application of T_EAP [m]
S = A_ref;                                                                  % Cross-section for N [m^2]

% Thrust imbalance

N = 1000;                                                                   % # samples
dT_max = 0.025;                                                             % Maximum thrust imbalance [%]
dT_1 = -dT_max + (dT_max+dT_max).*rand(N,1);                                % SRM left variation of thrust (rectangular distribution)
dT_2 = -dT_max + (dT_max+dT_max).*rand(N,1);                                % SRM right variation of thrust (rectangular distribution)

% Now are considered two cases in which is studied the effect of a normal
% aerodynamic force for different angles of attack

% Case 1: maximum T_SRM

n = 5;
alpha_agt = [0 3 5 7 10]*pi/180;
delta_EPC = zeros(n,N);
delta_EAP_1 = zeros(n,N);
delta_EAP_2 = zeros(n,N);
T_EPC_v = zeros(n,N);
T_EAP_1_v = zeros(n,N);
T_EAP_2_v = zeros(n,N);
T = zeros(n,N);
T_ideal = zeros(n,N);
options = optimoptions(@fsolve,'display','off');

for k = 1 : n
    % Deflection of the EPC nozzle to counteract N
    delta_EPC(k,:) = -asin(q(1)*S/(T_EPC(1)*CG(1))*...
                     abs(c_N_alpha(alpha_agt(k),0))*(CP_CG(1))*alpha_agt(k));
    % Random imbalance compensation by the SRMs
    for j = 1 : N
        if dT_1(j) > dT_2(j)
            f = @(d) -T_EAP(1)*(1 + dT_1(j))*cos(d)*y_bar + ...
                     T_EAP(1)*(1 + dT_1(j))*sin(d)*CG(1) + ...
                     T_EAP(1)*(1 + dT_2(j))*y_bar;
            d_1 = fsolve(f,0,options);
            delta_EAP_1(k,j) = -d_1;
        elseif dT_1(j) < dT_2(j)
            f = @(d) -T_EAP(1)*(1 + dT_1(j))*y_bar + ...
                     T_EAP(1)*(1 + dT_2(j))*cos(d)*y_bar - ...
                     T_EAP(1)*(1 + dT_2(j))*sin(d)*CG(1);
            d_1 = fsolve(f,0,options);
            delta_EAP_2(k,j) = d_1;
        end
        % Additional compensation of N with the two SRMs
        if delta_EPC(k,j) < -7*pi/180
            delta_EPC(k,j) = -7*pi/180;
            f = @(d) T_EPC(1)*sin(delta_EPC(1,j))*CG(2) + ...
                     T_EAP(1)*(1 + dT_1(j))*sin(d + delta_EAP_1(k,j))*CG(1) - ...
                     T_EAP(1)*(1 + dT_1(j))*cos(d + delta_EAP_1(k,j))*y_bar + ...
                     T_EAP(1)*(1 + dT_2(j))*sin(d + delta_EAP_2(k,j))*CG(1) + ...
                     T_EAP(1)*(1 + dT_2(j))*cos(d + delta_EAP_2(k,j))*y_bar - ...
                     q(1)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(1)*alpha_agt(k);
            d_1 = fsolve(f,0,options);
            delta_EAP_1(k,j) = delta_EAP_1(k,j) - d_1;
            delta_EAP_2(k,j) = delta_EAP_2(k,j) - d_1;
            f = @(d) T_EPC(1)*sin(delta_EPC(k,j))*CG(1) + ...
                            2*T_EAP(1)*sin(d)*CG(1) - ...
                            q(1)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(1)*alpha_agt(k);
                d_2 = fsolve(f,0,options);
        else
            d_2 = 0;
        end
        T_EPC_v(k,j) = T_EPC(1)*cos(delta_EPC(k,j));
        T_EAP_1_v(k,j) = T_EAP(1)*(1 + dT_1(j))*cos(delta_EAP_1(k,j));
        T_EAP_2_v(k,j) = T_EAP(1)*(1 + dT_2(j))*cos(delta_EAP_2(k,j));
        T_ideal(k,j) = T_EPC_v(k,j) + 2*T_EAP(1)*cos(d_2);
        T(k,j) = T_EPC_v(k,j) + T_EAP_1_v(k,j) + T_EAP_2_v(k,j);
    end
end

[dT,indexes] = sort(dT_1 - dT_2);
d_old = delta_EPC;
d_old_1 = delta_EAP_1;
d_old_2 = delta_EAP_2;
T_old = T;

for j = 1 : 1000
    delta_EPC(:,j) = d_old(:,indexes(j));
    delta_EAP_1(:,j) = d_old_1(:,indexes(j));
    delta_EAP_2(:,j) = d_old_2(:,indexes(j));
    T(:,j) = T_old(:,indexes(j));
end

figure()
plot(dT*100,delta_EPC(1,:)*180/pi,'b',...
     dT*100,delta_EPC(2,:)*180/pi,'r',...
     dT*100,delta_EPC(3,:)*180/pi,'g',...
     dT*100,delta_EPC(4,:)*180/pi,'c',...
     dT*100,delta_EPC(5,:)*180/pi,'y','linewidth',2.5)
grid on
title('Excursion at maximum T_E_A_P of the EPC nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_P_C [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','location','northeast')

figure()
hold on
plot(dT*100,delta_EAP_1(1,:)*180/pi,'ob',...
     dT*100,delta_EAP_1(2,:)*180/pi,'^r',...
     dT*100,delta_EAP_1(3,:)*180/pi,'sg',...
     dT*100,delta_EAP_1(4,:)*180/pi,'<c',...
     dT*100,delta_EAP_1(5,:)*180/pi,'>y')
hold off
grid on
title('Excursion at maximum T_E_A_P of the EAP left nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_A_P_1 [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','location','southwest')

figure()
hold on
plot(dT*100,delta_EAP_2(1,:)*180/pi,'ob',...
     dT*100,delta_EAP_2(2,:)*180/pi,'^r',...
     dT*100,delta_EAP_2(3,:)*180/pi,'sg',...
     dT*100,delta_EAP_2(4,:)*180/pi,'<c',...
     dT*100,delta_EAP_2(5,:)*180/pi,'>y')
hold off
grid on
title('Excursion at maximum T_E_A_P of the EAP right nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_A_P_2 [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','location','southwest')

figure()
hold on
plot(dT*100,T(1,:)/1e6,'ob',dT*100,T(2,:)/1e6,'^r',...
     dT*100,T(3,:)/1e6,'sg',dT*100,T(4,:)/1e6,'<c',...
     dT*100,T(5,:)/1e6,'>y')
plot(dT*100,T_ideal(1,:)/1e6,'b',...
     dT*100,T_ideal(2,:)/1e6,'r',...
     dT*100,T_ideal(3,:)/1e6,'g',...
     dT*100,T_ideal(4,:)/1e6,'c',...
     dT*100,T_ideal(5,:)/1e6,'y','linewidth',2.5)
hold off
grid on
title('Overall thrust at maximum T_E_A_P')
xlabel('\DeltaT [%]')
ylabel('T_t_o_t [MN]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
       '\alpha = 10°','T_i_d_e_a_l_1','T_i_d_e_a_l_2','T_i_d_e_a_l_3',...
       'T_i_d_e_a_l_4','T_i_d_e_a_l_5','location','northeast')

% Case 2: maximum q

n = 5;
alpha_agt = [0 3 5 7 10]*pi/180;
delta_EPC = zeros(n,N);
delta_EAP_1 = zeros(n,N);
delta_EAP_2 = zeros(n,N);
T_EPC_v = zeros(n,N);
T_EAP_1_v = zeros(n,N);
T_EAP_2_v = zeros(n,N);
T = zeros(n,N);
T_ideal = zeros(n,N);

for k = 1 : n
    % Deflection of the EPC nozzle to counteract N
    delta_EPC(k,:) = -asin(q(2)*S/(T_EPC(2)*CG(2))*...
                     abs(c_N_alpha(alpha_agt(k),0))*(CP_CG(2))*alpha_agt(k));
    % Random imbalance compensation by the SRMs
    for j = 1 : N
        if dT_1(j) > dT_2(j)
            f = @(d) -T_EAP(2)*(1 + dT_1(j))*cos(d)*y_bar + ...
                     T_EAP(2)*(1 + dT_1(j))*sin(d)*CG(2) + ...
                     T_EAP(2)*(1 + dT_2(j))*y_bar;
            d_1 = fsolve(f,0,options);
            delta_EAP_1(k,j) = -d_1;
        elseif dT_1(j) < dT_2(j)
            f = @(d) -T_EAP(2)*(1 + dT_1(j))*y_bar + ...
                     T_EAP(2)*(1 + dT_2(j))*cos(d)*y_bar - ...
                     T_EAP(2)*(1 + dT_2(j))*sin(d)*CG(2);
            d_1 = fsolve(f,0,options);
            delta_EAP_2(k,j) = d_1;
        end
        % Additional compensation of N with the two SRMs
        if delta_EPC(k,j) < -7*pi/180
            delta_EPC(k,j) = -7*pi/180;
            f = @(d) T_EPC(2)*sin(delta_EPC(1,j))*CG(2) + ...
                     T_EAP(2)*(1 + dT_1(j))*sin(d + delta_EAP_1(k,j))*CG(2) - ...
                     T_EAP(2)*(1 + dT_1(j))*cos(d + delta_EAP_1(k,j))*y_bar + ...
                     T_EAP(2)*(1 + dT_2(j))*sin(d + delta_EAP_2(k,j))*CG(2) + ...
                     T_EAP(2)*(1 + dT_2(j))*cos(d + delta_EAP_2(k,j))*y_bar - ...
                     q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
            d_1 = fsolve(f,0,options);
            delta_EAP_1(k,j) = delta_EAP_1(k,j) - d_1;
            delta_EAP_2(k,j) = delta_EAP_2(k,j) - d_1;
            f = @(d) T_EPC(2)*sin(delta_EPC(k,j))*CG(2) + ...
                            2*T_EAP(2)*sin(d)*CG(2) - ...
                            q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
                d_2 = fsolve(f,0,options);
        else
            d_2 = 0;
        end
        T_EPC_v(k,j) = T_EPC(2)*cos(delta_EPC(k,j));
        T_EAP_1_v(k,j) = T_EAP(2)*(1 + dT_1(j))*cos(delta_EAP_1(k,j));
        T_EAP_2_v(k,j) = T_EAP(2)*(1 + dT_2(j))*cos(delta_EAP_2(k,j));
        T_ideal(k,j) = T_EPC_v(k,j) + 2*T_EAP(2)*cos(d_2);
        T(k,j) = T_EPC_v(k,j) + T_EAP_1_v(k,j) + T_EAP_2_v(k,j);
    end
end

d_old = delta_EPC;
d_old_1 = delta_EAP_1;
d_old_2 = delta_EAP_2;
T_old = T;

for j = 1 : 1000
    delta_EPC(:,j) = d_old(:,indexes(j));
    delta_EAP_1(:,j) = d_old_1(:,indexes(j));
    delta_EAP_2(:,j) = d_old_2(:,indexes(j));
    T(:,j) = T_old(:,indexes(j));
end

figure()
plot(dT*100,delta_EPC(1,:)*180/pi,'b',...
     dT*100,delta_EPC(2,:)*180/pi,'r',...
     dT*100,delta_EPC(3,:)*180/pi,'g',...
     dT*100,delta_EPC(4,:)*180/pi,'c',...
     dT*100,delta_EPC(5,:)*180/pi,'y','linewidth',2.5)
grid on
title('Excursion at maximum q of the EPC nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_P_C [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','location','northeast')

figure()
hold on
plot(dT*100,delta_EAP_1(1,:)*180/pi,'ob',...
     dT*100,delta_EAP_1(2,:)*180/pi,'^r',...
     dT*100,delta_EAP_1(3,:)*180/pi,'sg',...
     dT*100,delta_EAP_1(4,:)*180/pi,'<c',...
     dT*100,delta_EAP_1(5,:)*180/pi,'>y')
patch([-5 5 5 -5],[-7.3 -7.3 -12 -12],'r','FaceAlpha',.25,'EdgeColor','none') 
hold off
grid on
title('Excursion at maximum q of the EAP left nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_A_P_1 [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','\delta_E_A_P_1 < -7.3° !','location','northeast')

figure()
hold on
plot(dT*100,delta_EAP_2(1,:)*180/pi,'ob',...
     dT*100,delta_EAP_2(2,:)*180/pi,'^r',...
     dT*100,delta_EAP_2(3,:)*180/pi,'sg',...
     dT*100,delta_EAP_2(4,:)*180/pi,'<c',...
     dT*100,delta_EAP_2(5,:)*180/pi,'>y')
patch([-5 5 5 -5],[-7.3 -7.3 -12 -12],'r','FaceAlpha',.25,'EdgeColor','none') 
hold off
grid on
title('Excursion at maximum q of the EAP right nozzle')
xlabel('\DeltaT [%]')
ylabel('\delta_E_A_P_2 [°]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
    '\alpha = 10°','\delta_E_A_P_2 < -7.3° !','location','northeast')
figure()
hold on
plot(dT*100,T(1,:)/1e6,'ob',dT*100,T(2,:)/1e6,'^r',...
     dT*100,T(3,:)/1e6,'sg',dT*100,T(4,:)/1e6,'<c',...
     dT*100,T(5,:)/1e6,'>y')
plot(dT*100,T_ideal(1,:)/1e6,'b',...
     dT*100,T_ideal(2,:)/1e6,'r',...
     dT*100,T_ideal(3,:)/1e6,'g',...
     dT*100,T_ideal(4,:)/1e6,'c',...
     dT*100,T_ideal(5,:)/1e6,'y','linewidth',2.5)
hold off
grid on
title('Overall thrust at maximum q')
xlabel('\DeltaT [%]')
ylabel('T_t_o_t [MN]')
legend('\alpha = 0°','\alpha = 3°','\alpha = 5°','\alpha = 7°',...
       '\alpha = 10°','T_i_d_e_a_l_1','T_i_d_e_a_l_2','T_i_d_e_a_l_3',...
       'T_i_d_e_a_l_4','T_i_d_e_a_l_5','location','northeast')