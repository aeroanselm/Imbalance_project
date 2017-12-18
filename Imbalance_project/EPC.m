% EPC.m - Ariane 5 ECA ---Imbalance project---
%  
% TYPE:
%   EPC auxiliary script
%
% DESCRIPTION: 
%   This routine is used to calculate the thrust profile of the EPC main core
%   stage.
%
% Notes for upgrading this script: 
%   Please do not upgrade directly this script. Make a your local copy and do
%   your modifications and upgrades and send me the upgraded script.
%   Contact the author for modifications.
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
%				MARCELLO SCALERA
%				ALESSIO NEGRI
%               NICOLÒ BRAMANI
%               ALBERTO SEDINI
%               PIERFRANCESCO PANSINI
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

[EPCout]=rocket(EPCdata);

EPC_Pe     = EPCout.Pe;    
EPC_Ve     = EPCout.Ve;    
EPC_rho_cc = EPCout.rho_cc;
EPC_rho_e  = EPCout.rho_e;
EPC_mdot_p = EPCout.mdot_p;
EPC_A_e    = EPCout.A_e;   
EPC_a_e    = EPCout.a_e;   
EPC_Te     = EPCout.Te;    
EPC_Me     = EPCout.Me;    
EPC_IspSL  = EPCout.IspSL; 
EPC_IspVac = EPCout.IspVac;

a = 0.0065;
g = 9.81;
T_z0 = 288.16;
T_z11000 = T_z0-a*11000;
P_z0 = 101325;
R = 287;
n = g/(R*a);

P_amb_z11000 = P_z0*(1-a*11000/T_z0).^n;


z = [0:145000];
P_amb = [P_z0*(1-a*z(1:11001)/T_z0).^n P_amb_z11000*exp(-g/(R*T_z11000)*(z(11002:145001)-11000))];

EPC_P_amb = P_amb(1:145001);

EPC_T = (EPC_mdot_p*EPC_Ve + (EPC_Pe-EPC_P_amb)*EPC_A_e);

figure()
hold on
grid on
plot(z(1:145001),EPC_T)
xlabel('Altitude [km]')
ylabel('Thrust [N]')
title('EPC Thrust vs Altitude')

