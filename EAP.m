% EAP.m - Ariane 5 ECA ---Imbalance project---
%  
% TYPE:
%   EAP auxiliary script
%
% DESCRIPTION: 
%   This routine is used to calculate the thrust profile of the EAP solid
%   rocket boosters. In order to do this, linear interpolation method is used,
%   starting from graphical data computed using the document A5prop.pdf, provided
%   along with this script.
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

EAP_T = @(t) interp1(x,y,t);
EAP_mdot = @(t) EAP_T(t)/g0/EAP_Isp;

figure()
hold on
grid on
plot(0:130,EAP_T(0:130));           % for debugging only
xlabel('Time [sec]')
ylabel('Thrust [N]')
title('EAP Thrust vs Time')
