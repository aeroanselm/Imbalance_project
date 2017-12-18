close all
clear
clc

% Aerodynamics and Thrust Vector Control (TVC) of Ariane 5-ECA

%% AERODYNAMICS

% DATA

l_B = 50.5;                                                                 % Body length [m]
l_N = 7;                                                                    % Nose length [m]
l_F = 17;                                                                   % Fairing length [m]
l_ESCA = 4.711;                                                             % ESC-A length [m]
l_EPC = 23.8;                                                               % EPC length [m]
l_SRM = 30;                                                                 % SRM length [m]
l_N_SRM = 4.375;                                                            % SRM nose length [m]
x_f = 28.98;                                                                % SRM flare position [m]
d_c = 1;                                                                    % Circle length [m]
d_ESCA = 5.4;                                                               % ESC-A diameter [m]
d_EPC = d_ESCA;                                                             % EPC diameter [m]
d_SRM = 3.05;                                                               % SRM diameter [m]
d_f = 3.10;                                                                 % SRM flare diameter [m]
h_f = 1.02;                                                                 % SRM flare heigth [m]
d_m = 0.5*(d_SRM+d_f);                                                      % SRM flare mean diameter [m]
M_F = 2675 + 10000;                                                         % Fairing + Payload mass [kg]
M_ESCA = 4540 + 14900;                                                      % ESC-A mass [kg]
M_EPC = 12500 + 133000 + 26000;                                             % EPC mass [kg]
M_SRM = 37000+238000;                                                       % SRM mass [kg]
M_tot = M_F+M_ESCA+M_EPC+2*M_SRM;                                           % Ariane 5-ECA mass [kg]

% CALCULATIONS of areas

a_2 = d_EPC+2*d_SRM;                                                        % Major asis [m]
a = a_2/2;                                                                  % Semi-major axis [m]
b = d_ESCA/2;                                                               % Semi-minor axis [m]
A_nose = pi/4*d_c^2;                                                        % Area of the fictitious hemispherical nose [m^2]
A_SRM = pi/4*d_SRM^2;                                                       % SRM area [m^2]
A_EPC = pi/4*d_EPC^2;                                                       % EPC area [m^2]
A_ref = A_EPC + 2*A_SRM;                                                    % Reference area [m^2]
d_e_Vulcain = 2.10;                                                         % Vulcain exit nozzle diameter [m]
A_e_EPC = pi/4*d_e_Vulcain^2;                                               % Vulcain exit area [m^2]
d_e_SRM = 3.10;                                                             % SRM exit nozzle diameter [m]
A_e_SRM = pi/4*d_e_SRM^2;                                                   % SRM exit area [m^2]
A_e_bar = A_e_EPC +2*A_e_SRM;                                               % Sum of all the exit nozzle areas [m^2]
A_b = pi*a*b;                                                               % Elliptical base area [m^2]
A_f = pi/4*d_f^2;                                                           % Flare base area [m^2]

% C_D CALCULATIONS

%load('cd_Ma.m')
%c_D_data = c_D_vs_Ma(:,2);                                                  % c_D real data
%Ma_data = c_D_vs_Ma(:,1);                                                   % Ma real data
Ma = linspace(0,6,1000);
%c_D = spline(Ma_data,c_D_data,Ma);

figure(1)
%plot(Ma,cd_Ma(Ma),'b')
grid on
ylim([0,1])
title('Ariane 5-ECA drag coefficient')
xlabel('M')
ylabel('c_D')

% C_N & C_N_ALPHA CALCULATIONS

alpha_v = linspace(0,pi/2,1000);
c_N = @(alpha,phi) (a/b*cos(phi)^2 +...
    b/a*sin(phi)^2)*(abs(sin(2*alpha).*cos(alpha/2))+ ...
    1.3*l_B/(2*sqrt(a*b))*sin(alpha).^2);
c_N_alpha = @(alpha,phi) (a/b*cos(phi)^2 +...
    b/a*sin(phi)^2)*(abs(2*cos(2*alpha).*cos(alpha/2)-...
    0.5*sin(2*alpha).*sin(alpha./2)) +...
    1.3*l_B/(2*sqrt(a*b))*2*sin(alpha).*cos(alpha));

figure(2)
plot(alpha_v*180/pi,c_N(alpha_v,0),'b',...
     alpha_v*180/pi,c_N(alpha_v,5*pi/180),'r',...
     alpha_v*180/pi,c_N(alpha_v,10*pi/180),'g')
grid on
title('Normal force coefficient')
xlabel('\alpha [�]')
ylabel('c_N')
legend('\phi = 0�','\phi = 5�','\phi = 10�','location','northwest')

figure(3)
plot(alpha_v*180/pi,c_N_alpha(alpha_v,0),'b',...
     alpha_v*180/pi,c_N_alpha(alpha_v,5*pi/180),'r',...
     alpha_v*180/pi,c_N_alpha(alpha_v,10*pi/180),'g')
grid on
title('Normal force coefficient derived by \alpha')
xlabel('\alpha [�]')
ylabel('c_N_/_\alpha')
legend('\phi = 0�','\phi = 5�','\phi = 10�','location','northwest')

% X_CP CALCULATIONS

x_CP_mb = 2/3*l_N;
x_CP_SRM = 2/3*l_N_SRM*A_SRM/A_f + (1-A_SRM/A_f)*l_SRM - h_f*((d_m/d_SRM)^2 - 1)*A_SRM/A_f;
x_CP = (x_CP_mb*(l_B-l_SRM) + (x_CP_SRM + l_B-l_SRM)*2*l_SRM)/(l_B-l_SRM+2*l_SRM);

%% TVC

% Data

CG = [18.1752 17.5384 17.0891];                                             % CG positions from the base [m]
CP_CG = l_B - x_CP - CG;                                                    % CP-CG positions from the base [m]
% 18 s - 71 s - 97 s
% 394.7 m - 14570 m  - 29030 m
T_SRM = [6485000 5242000];% 5997000];
T_EPC = [934500 1251000];% 1294000];
q = [4686 34000];% 12370];
y_bar = 2.7 + 1.525;
S = A_ref;

% Thrust imbalance

N = 1000;                                                                   % # samples
dT_max = 0.025;
dT_1 = -dT_max + (dT_max+dT_max).*rand(N,1);                                % SRM left variation of thrust
dT_2 = -dT_max + (dT_max+dT_max).*rand(N,1);                                % SRM right variation of thrust
N_bins = floor(1.87*(N-1)^0.4 + 1);
figure(4)
plot(linspace(1,N,N),dT_1,'ob',linspace(1,N,N),dT_2,'or')
grid on

% Case 1: vertical launch at maximum T_SRM

delta_SRM_1 = zeros(1,N);
delta_SRM_2 = zeros(1,N);
T_SRM_1_v = zeros(1,N);
T_SRM_2_v = zeros(1,N);
T = zeros(1,N);
T_ideal = (T_EPC(1)+2*T_SRM(1))*ones(N,1);
options = optimoptions(@fsolve,'display','off');
for j = 1 : N % Random imbalance
    if dT_1(j) > dT_2(j)
        f = @(d) -T_SRM(1)*(1 + dT_1(j))*cos(d)*y_bar + ...
                 T_SRM(1)*(1 + dT_1(j))*sin(d)*CG(1) + ...
                 T_SRM(1)*(1 + dT_2(j))*y_bar;
        d_1 = fsolve(f,0,options);
        delta_SRM_1(j) = -d_1;
    elseif dT_1(j) < dT_2(j)
        f = @(d) -T_SRM(1)*(1 + dT_1(j))*y_bar + ...
                 T_SRM(1)*(1 + dT_2(j))*cos(d)*y_bar - ...
                 T_SRM(1)*(1 + dT_2(j))*sin(d)*CG(1);
        d_1 = fsolve(f,0,options);
        delta_SRM_2(j) = d_1;
    end
    T_SRM_1_v(j) = T_SRM(1)*(1 + dT_1(j))*cos(delta_SRM_1(j));
    T_SRM_2_v(j) = T_SRM(1)*(1 + dT_2(j))*cos(delta_SRM_2(j));
    T(j) = T_EPC(1) + T_SRM_1_v(j) + T_SRM_2_v(j);
end
[dT,indexes] = sort(dT_1 - dT_2);
d_old_1 = delta_SRM_1;
d_old_2 = delta_SRM_2;
T_old = T;
for j = 1 : 1000
    delta_SRM_1(j) = d_old_1(indexes(j));
    delta_SRM_2(j) = d_old_2(indexes(j));
    T(j) = T_old(indexes(j));
end
figure(5)
subplot(2,1,1)
plot(dT*100,delta_SRM_1*180/pi,'or')
grid on
title('Escursion during vertical path of the SRM left')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_1 [�]')
subplot(2,1,2)
plot(dT*100,delta_SRM_2*180/pi,'og')
grid on
title('Escursion during vertical path of the SRM right')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_2 [�]')
figure(6)
plot(dT*100,T_ideal/1e6,'k',dT*100,T/1e6,'om')
grid on
title('Overall thrust during vertical path')
xlabel('\DeltaT [%]')
ylabel('T_t_o_t [MN]')
legend('T_i_d_e_a_l','location','northeast')

% Case 2: maximum q

n = 3;
alpha_agt = [1 5 7]*pi/180;
delta_EPC = zeros(n,N);
delta_SRM_1 = zeros(n,N);
delta_SRM_2 = zeros(n,N);
T_EPC_v = zeros(n,N);
T_SRM_1_v = zeros(n,N);
T_SRM_2_v = zeros(n,N);
T = zeros(n,N);
T_ideal = zeros(n,N);
for k = 1 : n
    delta_EPC(k,:) = -asin(q(2)*S/(T_EPC(2)*CG(2))*...
                     abs(c_N_alpha(alpha_agt(k),0))*(CP_CG(2))*alpha_agt(k));
    for j = 1 : N % Random imbalance
        if dT_1(j) > dT_2(j)
            f = @(d) -T_SRM(2)*(1 + dT_1(j))*cos(d)*y_bar + ...
                     T_SRM(2)*(1 + dT_1(j))*sin(d)*CG(2) + ...
                     T_SRM(2)*(1 + dT_2(j))*y_bar;
            d_1 = fsolve(f,0,options);
            delta_SRM_1(k,j) = -d_1;
        elseif dT_1(j) < dT_2(j)
            f = @(d) -T_SRM(2)*(1 + dT_1(j))*y_bar + ...
                     T_SRM(2)*(1 + dT_2(j))*cos(d)*y_bar - ...
                     T_SRM(2)*(1 + dT_2(j))*sin(d)*CG(2);
            d_1 = fsolve(f,0,options);
            delta_SRM_2(k,j) = d_1;
        end
        if delta_EPC(k,j) < -7*pi/180
            delta_EPC(k,j) = -7*pi/180;
            f = @(d) T_EPC(2)*sin(delta_EPC(1,j))*CG(2) + ...
                     T_SRM(2)*(1 + dT_1(j))*sin(d + delta_SRM_1(k,j))*CG(2) - ...
                     T_SRM(2)*(1 + dT_1(j))*cos(d + delta_SRM_1(k,j))*y_bar + ...
                     T_SRM(2)*(1 + dT_2(j))*sin(d + delta_SRM_2(k,j))*CG(2) + ...
                     T_SRM(2)*(1 + dT_2(j))*cos(d + delta_SRM_2(k,j))*y_bar - ...
                     q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
            d_1 = fsolve(f,0,options);
            delta_SRM_1(k,j) = delta_SRM_1(k,j) - d_1;
            delta_SRM_2(k,j) = delta_SRM_2(k,j) - d_1;
            f = @(d) T_EPC(2)*sin(delta_EPC(k,j))*CG(2) + ...
                            2*T_SRM(2)*sin(d)*CG(2) - ...
                            q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
                d_2 = fsolve(f,0,options);
        else
            d_2 = 0;
        end
        T_EPC_v(k,j) = T_EPC(2)*cos(delta_EPC(k,j));
        T_SRM_1_v(k,j) = T_SRM(2)*(1 + dT_1(j))*cos(delta_SRM_1(k,j));
        T_SRM_2_v(k,j) = T_SRM(2)*(1 + dT_2(j))*cos(delta_SRM_2(k,j));
        T_ideal(k,j) = T_EPC_v(k,j) + 2*T_SRM(2)*cos(d_2);
        T(k,j) = T_EPC_v(k,j) + T_SRM_1_v(k,j) + T_SRM_2_v(k,j);
    end
end
d_old = delta_EPC;
d_old_1 = delta_SRM_1;
d_old_2 = delta_SRM_2;
T_old = T;
for j = 1 : 1000
    delta_EPC(:,j) = d_old(:,indexes(j));
    delta_SRM_1(:,j) = d_old_1(:,indexes(j));
    delta_SRM_2(:,j) = d_old_2(:,indexes(j));
    T(:,j) = T_old(:,indexes(j));
end
%%
figure(7)
plot(dT*100,delta_EPC(1,:)*180/pi,'ob',...
     dT*100,delta_EPC(2,:)*180/pi,'^r',...
     dT*100,delta_EPC(3,:)*180/pi,'sg')
grid on
title('Escursion at maximum q of the EPC')
xlabel('\DeltaT [%]')
ylabel('\delta_E_P_C [�]')
legend('\alpha = 1�','\alpha = 3�','\alpha = 5�',...
    'location','northeast')
ylim([-7.5,-3.5])
figure(8)
hold on
plot(dT*100,delta_SRM_1(1,:)*180/pi,'ob',...
     dT*100,delta_SRM_1(2,:)*180/pi,'^r',...
     dT*100,delta_SRM_1(3,:)*180/pi,'sg')
patch([-5 5 5 -5],[-7.3 -7.3 -9 -9],'r','FaceAlpha',.5,'EdgeColor','none') 
hold off
grid on
title('Escursion at maximum q of the SRM left')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_1 [�]')
legend('\alpha = 1�','\alpha = 3�','\alpha = 5�',...
    '\delta_S_R_M_1 < -7.3� !','location','northwest')
figure(9)
plot(dT*100,delta_SRM_2(1,:)*180/pi,'ob',...
     dT*100,delta_SRM_2(2,:)*180/pi,'^r',...
     dT*100,delta_SRM_2(3,:)*180/pi,'sg')
grid on
title('Escursion at maximum q of the SRM right')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_2 [�]')
legend('\alpha = 1�','\alpha = 3�','\alpha = 5�',...
    'location','southeast')
figure(10)
hold on
plot(dT*100,T(1,:)/1e6,'ob',dT*100,T(2,:)/1e6,'^r',...
    dT*100,T(3,:)/1e6,'sg')
plot(dT*100,T_ideal(1,:)/1e6,'b','linewidth',2.5)
plot(dT*100,T_ideal(2,:)/1e6,'r','linewidth',2.5)
plot(dT*100,T_ideal(3,:)/1e6,'g','linewidth',2.5)
hold off
grid on
title('Overall thrust at maximum q')
xlabel('\DeltaT [%]')
ylabel('T_t_o_t [MN]')
legend('T_i_d_e_a_l','location','northeast')
legend('\alpha = 1�','\alpha = 3�','\alpha = 5�',...
    'location','northeast')
% figure(7)
% plot(dT*100,delta_EPC*180/pi,'ob')
% grid on
% title('Escursion at maximum q of the EPC')
% xlabel('\DeltaT [%]')
% ylabel('\delta_E_P_C [�]')
% figure(8)
% subplot(2,1,1)
% plot(dT*100,delta_SRM_1*180/pi,'or')
% grid on
% title('Escursion at maximum q of the SRM left')
% xlabel('\DeltaT [%]')
% ylabel('\delta_S_R_M_1 [�]')
% subplot(2,1,2)
% plot(dT*100,delta_SRM_2*180/pi,'og')
% grid on
% title('Escursion at maximum q of the SRM right')
% xlabel('\DeltaT [%]')
% ylabel('\delta_S_R_M_2 [�]')
% figure(9)
% plot(dT*100,T_ideal/1e6,'k',dT*100,T/1e6,'om')
% grid on
% title('Overall thrust at maximum q')
% xlabel('\DeltaT [%]')
% ylabel('T_t_o_t [MN]')
% legend('T_i_d_e_a_l','location','northeast')

% figure(10)
% plot(T_SRM(3)*dT/1e6,delta_EPC(2,:)*180/pi,'b')
% grid on
% title('Escursion at maximum q of the LRE')
% xlabel('\DeltaT [MN]')
% ylabel('\delta_E_P_C [�]')
% figure(11)
% subplot(2,1,1)
% plot(T_SRM(3)*dT/1e6,delta_SRM_1(2,:)*180/pi,'r')
% grid on
% title('Escursion at maximum q of the SRM left')
% xlabel('\DeltaT [MN]')
% ylabel('\delta_S_R_M_1 [�]')
% subplot(2,1,2)
% plot(T_SRM(3)*dT/1e6,delta_SRM_2(2,:)*180/pi,'g')
% grid on
% title('Escursion at maximum q of the SRM right')
% xlabel('\DeltaT [MN]')
% ylabel('\delta_S_R_M_2 [�]')
% figure(12)
% plot(T_SRM(3)*dT/1e6,T_ideal(2,:)/1e6,'k',T_SRM(3)*dT/1e6,T(2,:)/1e6,'m')
% grid on
% title('Overall thrust at maximum q')
% xlabel('\DeltaT [MN]')
% ylabel('T_t_o_t [MN]')
% legend('T_i_d_e_a_l','location','northeast')
return
% Case 3: maximum q and random alpha

n = 3;
delta_EPC = zeros(n,N);
delta_SRM_1 = zeros(n,N);
delta_SRM_2 = zeros(n,N);
T_EPC_v = zeros(n,N);
T_SRM_1_v = zeros(n,N);
T_SRM_2_v = zeros(n,N);
T = zeros(n,N);
T_ideal = zeros(n,N);
alpha_agt = -[10 12.5 15]*pi/180;
for k = 1 : n % Changing alpha
    delta_EPC(k,:) = -q(2)*S/(T_EPC(2)*CG(2))*...
        abs(c_N_alpha(alpha_agt(k),0))*(CP_CG(2))*alpha_agt(k);
    for j = 1 : N % Random imbalance
        if delta_EPC(k,j) > 7.3*pi/180
            delta_EPC(k,j) = 7.3*pi/180;
            f = @(d) T_EPC(2)*sin(delta_EPC(k,j))*CG(2) + ...
                        2*T_SRM(2)*sin(d)*CG(2) + ...
                        q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
            d_1 = fsolve(f,0,options);
            delta_SRM_1(k,j) = d_1;
            delta_SRM_2(k,j) = d_1;
            d_2 = d_1;
        elseif delta_EPC(k,j) < -7.3*pi/180
            delta_EPC(k,j) = -7.3*pi/180;
            f = @(d) T_EPC(2)*sin(delta_EPC(k,j))*CG(2) + ...
                        2*T_SRM(2)*sin(d)*CG(2) + ...
                        q(2)*S*abs(c_N_alpha(alpha_agt(k),0))*CP_CG(2)*alpha_agt(k);
            d_1 = fsolve(f,0,options);
            delta_SRM_1(k,j) = d_1;
            delta_SRM_2(k,j) = d_1;
            d_2 = d_1;
        else
            d_2 = 0;
        end
        f = @(d) -T_SRM(2)*(1 + dT(j))*cos(d)*y_bar + ...
            T_SRM(2)*(1 + dT(j))*sin(d)*CG(2) + ...
            T_SRM(2)*y_bar;
        d_1 = fsolve(f,0,options);
        if dT(j) >= 0
            delta_SRM_1(k,j) = delta_SRM_1(k,j) + d_1;
        else
            delta_SRM_2(k,j) = delta_SRM_2(k,j) + d_1;
        end
        T_EPC_v(k,j) = T_EPC(2)*cos(delta_EPC(k,j));
        T_SRM_1_v(k,j) = T_SRM(2)*cos(delta_SRM_1(k,j));
        T_SRM_2_v(k,j) = T_SRM(2)*cos(delta_SRM_2(k,j));
        T_ideal(k,j) = T_EPC_v(k,j) + 2*T_SRM(2)*cos(d_2);
        T(k,j) = T_EPC_v(k,j) + T_SRM_1_v(k,j) + T_SRM_2_v(k,j);
    end
end
%%
figure(13)
plot(dT*100,delta_EPC(1,:)*180/pi,'ob',...
     dT*100,delta_EPC(2,:)*180/pi,'^r',...
     dT*100,delta_EPC(3,:)*180/pi,'sg')
grid on
title('Escursion at maximum q of the EPC')
xlabel('\DeltaT [%]')
ylabel('\delta_E_P_C [�]')
legend('\alpha = - 10�','\alpha = - 12.5�','\alpha = - 15�',...
    'location','northeast')
figure(14)
hold on
plot(dT*100,delta_SRM_1(1,:)*180/pi,'ob',...
     dT*100,delta_SRM_1(2,:)*180/pi,'^r',...
     dT*100,delta_SRM_1(3,:)*180/pi,'sg')
patch([-15 15 15 -15],[7.3 7.3 9 9],'r','FaceAlpha',.5,'EdgeColor','none') 
hold off
grid on
title('Escursion at maximum q of the SRM left')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_1 [�]')
legend('\alpha = - 10�','\alpha = - 12.5�','\alpha = - 15�',...
    '\delta_S_R_M_1 > 7.3� !','location','northwest')
figure(15)
plot(dT*100,delta_SRM_2(1,:)*180/pi,'ob',...
     dT*100,delta_SRM_2(2,:)*180/pi,'^r',...
     dT*100,delta_SRM_2(3,:)*180/pi,'sg')
grid on
title('Escursion at maximum q of the SRM right')
xlabel('\DeltaT [%]')
ylabel('\delta_S_R_M_2 [�]')
legend('\alpha = - 10�','\alpha = - 12.5�','\alpha = - 15�',...
    'location','southeast')
figure(16)
plot(dT*100,T(1,:)/1e6,'ob',dT*100,T(2,:)/1e6,'^r',...
    dT*100,T(3,:)/1e6,'sg',dT*100,T_ideal(1,:)/1e6,'k',...
    dT*100,T_ideal(2,:)/1e6,'k',...
    dT*100,T_ideal(3,:)/1e6,'k')
grid on
title('Overall thrust at maximum q')
xlabel('\DeltaT [%]')
ylabel('T_t_o_t [MN]')
legend('T_i_d_e_a_l','location','northeast')
legend('\alpha = - 10�','\alpha = - 12.5�','\alpha = - 15�',...
    'location','southeast')