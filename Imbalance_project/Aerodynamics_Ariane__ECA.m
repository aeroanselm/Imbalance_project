% CALCULATIONS

d_b = d_EPC+2*d_EAP;                                                        % Base diameter [m]
a = d_b/2;
b = d_ESCA/2;
A_ref = pi/4*d_ref^2;
A_nose = pi/4*d_c^2;
A_EAP = pi/4*d_EAP^2;
d_e_Vulcain = 2.10;
A_e_EPC = pi/4*d_e_Vulcain^2;
d_e_EAP = 3.10;
A_e_EAP = pi/4*d_e_EAP^2;
A_e_bar = A_e_EPC +2*A_e_EAP;
A_ref_bar = A_ref + 2*A_EAP;
A_b = pi/4*d_b^2;
x_F = l_F/2;
x_ESCA = l_F + l_ESCA/2;
x_EPC = l_F + l_ESCA + l_EPC/2;
x_EAP = (l_B-l_EAP) + l_EAP/2;


% C_D CALCULATIONS

t = linspace(0,141.9,1000);
t_data = [0 12.77 17.05 22.60 32.05 48.96 68.64 141.9]; % time [s]
h_data = [0 90 320 870 2410 6690 13500 66700]; % altitude [m]
v_data = [0 36.2 72.2 123.1 209.9 322.3 523.6 2020.2]; % velocity [m/s]
h = spline(t_data,h_data,t);
v = spline(t_data,v_data,t);
figure(1)
subplot(2,2,1)
plot(t,h,'r',t_data,h_data,'o')
grid on
title('Altitude')
subplot(2,2,2)
plot(t,v,'r',t_data,v_data,'o')
grid on
title('Velocity')
rho = zeros(1,1000);
for k = 1 : 1000
    rho(k) = EDMA(h(k)/1000);
end
subplot(2,2,3)
plot(t,rho,'r')
grid on
title('Density')
q = 0.5.*rho.*(v.^2);
subplot(2,2,4)
plot(t,q,'r')
grid on
title('Dynamic ressure')
% Ma = 1 at 49 s , q_max at 1 min 08 s , burnout at 2 min 21 s
% q = 35000/0.02216*normpdf(linspace(0,180,1000),68,18);
% vect = linspace(-5,1,1000);
% v = 2000/exp(1)*exp(vect);
Ma = linspace(0,5,1000);
% vect_Ma = linspace(1,10,1000);
% Ma = log(vect_Ma);
% rho = 2*q./(v.^2);

c_D_0_bw = zeros(1,1000);
for k = 1 : length(Ma)
    if Ma(k) <= 1
        c_D_0_bw(k) = 0;
%     elseif Ma(k) >= 5
%         c_D_0_bw(k) = 1.586*(atan(0.5/l_N/d_ESCA))^1.69*...
%             (A_ref-A_nose)/A_ref + 0.665*1.586*A_nose/A_ref;
    else
        c_D_0_bw(k) = (1.586+1.834/(Ma(k)^2))*(atan(0.5/l_N/d_ESCA))^1.69*...
            (A_ref-A_nose)/A_ref + ...
            0.665*(1.586+1.834/(Ma(k)^2))*A_nose/A_ref;
        c_D_0_bw(k) = c_D_0_bw(k) + 2*(1.586+1.834/(Ma(k)^2))*...
            (atan(0.5/l_N_EAP/d_EAP))^1.69*A_EAP/A_ref;
    end
end

c_D_0_bf = zeros(1,1000);
% h = linspace(0,100,1000);
% gamma = 1.4;
% P = 101325.*exp(-9.80665*0.0289644.*h/8.31447/288.15);
% rho = zeros(1,1000);
% for k = 1 : length(Ma)
%     rho(k) = atmosphere(h(k));
% end
vect = linspace(-5,1,1000);
% q = 0.5*rho.*v.^2;
% R = 8.31447/0.0289644;
% T = P./(R*rho);
% a_s = sqrt(gamma*R*T);
% q = 0.5*rho.*a_s.^2.*Ma.^2;
% h = 13.567;
% q = 37.5e3;

for k = 1 : length(Ma)
    c_D_0_bf(k) = 0.053*l_B/d_ESCA*(Ma(k)/(q(k)*0.0209)/(l_B*3.28))^0.2...
        + 2*0.053*l_EAP/d_EAP*(Ma(k)/(q(k)*0.0209)/(l_EAP*3.28))^0.2;
end
% c_D_0_bf = zeros(1,1000);
c_D_0_b = zeros(1,1000);
for k = 1 : length(Ma)
    if Ma(k)<= 1
        c_D_0_b(k) = (1-A_e_bar/A_ref_bar)*(0.12+0.13*Ma(k)^2)*...
            A_ref_bar/A_ref;
    else
        c_D_0_b(k) = (1-A_e_bar/A_ref_bar)*0.25/Ma(k)*...
            A_ref_bar/A_ref;
    end
end

c_D_0_body = c_D_0_bw + c_D_0_bf + c_D_0_b;

figure(2)
hold on
plot(Ma,c_D_0_bw,'--b')
plot(Ma,c_D_0_bf,'--r')
plot(Ma,c_D_0_b,'--g')
plot(Ma,c_D_0_body,'-.k')
hold off
grid on
title('Drag coefficient')
xlabel('Ma')
ylabel('c_D_0^b^o^d^y')
legend('c_D_0^b^o^d^y^ ^w^a^v^e','c_D_0^b^o^d^y^ ^f^r^i^c^t^i^o^n',...
    'c_D_0^b^a^s^e','c_D_0^b^o^d^y','location','northeast')

% C_N CALCULATIONS

alpha = linspace(0,pi/2,1000);
phi = 0*pi/180;
c_N = ((abs(sin(2*alpha).*cos(alpha/2))+1.3*l_B/d_ESCA*sin(alpha).^2)*l_B + ...
    (a/b*cos(phi)^2 +b/a*sin(phi)^2)*(abs(sin(2*alpha).*cos(alpha/2))+ ...
    1.3*l_EAP/(2*(a*b)^0.5)*sin(alpha).^2)*A_b/A_ref*l_EAP)/(l_B+l_EAP);
c_N_alpha = 2*A_b/A_ref;
figure(3)
plot(alpha*180/pi,c_N,'-b')
grid on
title('Body normal force')
xlabel('alpha [°]')
ylabel('c_N')

% E CALCULATIONS

c_D_0 = 0.2;
E = (c_N.*cos(alpha)-c_D_0.*sin(alpha))./...
    (c_N.*sin(alpha)+c_D_0.*cos(alpha));
figure(4)
plot(alpha*180/pi,E,'-b')
grid on
title('Aerodynamic efficiency')
xlabel('alpha [°]')
ylabel('E')

% STABILITY

x_CG = (x_F*M_F + x_ESCA*M_ESCA + x_EPC*M_EPC + 2*x_EAP*M_EAP)/...
    (M_F + M_ESCA + M_EPC + 2*M_EAP);
fprintf('Position of the center of gravity x_CG = %f [m]\n\n',x_CG)
x_CP = d_ref*(2/3*l_N/d_ESCA + (1-A_ref/A_b)*l_B/d_ESCA);
fprintf('Position of the center of pressure x_CP = %f [m]\n\n',x_CP)
m_s = (x_CG-x_CP)/d_ESCA;
fprintf('Margin of stability %f\n\n',m_s)
if m_s < 0
    disp('The launcher is statically stable')
    fprintf('\n')
else
    disp('The launcher is statically unstable')
    fprintf('\n')
end