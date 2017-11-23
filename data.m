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

