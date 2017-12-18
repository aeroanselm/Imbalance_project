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


% Aerodynamics

