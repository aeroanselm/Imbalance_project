% data.m - Ariane 5 ECA ---Imbalance project---
%  
% TYPE:
%   Data Script
%
% DESCRIPTION: 
%   This routine is used to import data relative to the Ariane 5 ECA 
%   launcher by ESA.  
%       
% Notes for upgrading this script: 
%   Do not change the structure of the script. In this script structures are
%   used in order to provide in the correct way data to the functions used 
%   in other scripts of this project. DO NOT add constants that can be easily
%   computed starting form other ones (avoid redundancy).
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

% DEFINED CONSTANTS
g0 = 9.81;

% CRYOGENIC MAIN CORE STAGE (EPC)
EPCdata = struct('mprop',   173300,  ...
                 'e_ratio', 61.5,    ...
                 'gamma',   1.2873,  ...
                 'Pcc',    11600000,...
                 'Tcc',    3539.57, ...
                 'mmol',    13.534,  ...
                 'mp',      175000,  ...
                 'tburn',   540,     ...
                 'nengines',1);

% SOLID ROCKET BOOSTER (EAP)
    
load('thrust_b.mat');
t = 0:130;
EAP_IspSL = 262;
EAP_IspVac = 274.5;
EAP_Isp = (EAP_IspSL+EAP_IspVac)/2;

% GEOMETRICAL AND STRUCTURAL DATA
l_B = 50.5;                                                                 % Body length [m]
l_N = 7;                                                                    % Nose length [m]
l_F = 17;                                                                   % Fairing length [m]
l_ESCA = 4.711;                                                             % ESC-A length [m]
l_EPC = 23.8;                                                               % EPC length [m]
l_EAP = 30;                                                                 % EAP length [m]
l_N_EAP = 4.375;                                                            % EAP nose length [m]
x_f = 28.98;                                                                % EAP flare position [m]
d_c = 1;                                                                    % Circle length [m]
d_ESCA = 5.4;                                                               % ESC-A diameter [m]
d_EPC = d_ESCA;                                                             % EPC diameter [m]
d_EAP = 3.05;                                                               % EAP diameter [m]
d_f = 3.10;                                                                 % EAP flare diameter [m]
h_f = 1.02;                                                                 % EAP flare heigth [m]
d_m = 0.5*(d_EAP+d_f);                                                      % EAP flare mean diameter [m]
M_F = 2675 + 10000;                                                         % Fairing + Payload mass [kg]
M_ESCA = 4540 + 14900;                                                      % ESC-A mass [kg]
M_EPC = 12500 + 133000 + 26000;                                             % EPC mass [kg]
M_EAP = 37000+238000;                                                       % EAP mass [kg]
M_tot = M_F+M_ESCA+M_EPC+2*M_EAP;                                           % Ariane 5-ECA mass [kg]

