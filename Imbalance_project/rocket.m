% TYPE:
%	[function]
% PROTOTYPE:
%	[out] = rocket(data)
% FUTURE DEVELOPMENT:
%   Work in progress.
% ORIGINAL VERSION:
%   08/01/2017, ALESSANDRO MARIA MASSERINI
% ------------------------------------------------------
% ***** BEGIN GPL LICENSE BLOCK *****
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the improplied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ***** END GPL LICENCE BLOCK *****
% ------------------------------------------------------
% AUTHORS:
%   Name: ALESSANDRO MARIA 
%   Surname: MASSERINI
%   ID number: 808451, 883251
%   Course: Aerospace Engineering
%           Space Engineering
%   Department: DAER
%   University: Politecnico di Milano
%   Class: Aerospace Propulsion
%   Creation: 08/01/2017
%   Contact: alessandro.masserini@mail.polimi.it
% CHANGELOG: 
%   ALESSANDRO MARIA MASSERINI creation of this function - 08/01/2017
%   ALESSANDRO MARIA MASSERINI substituting inp and out variables with structs - 23/11/2017                           
%%------------------------------------------------------ 
function [out] = rocket(data)

e_ratio   = data.e_ratio  ;
gamma     = data.gamma    ;
Pcc       = data.Pcc      ;
Tcc       = data.Tcc      ;
mmol      = data.mmol     ;
mprop     = data.mprop    ;
tburn     = data.tburn    ;
nengines  = data.nengines ;

R = 8314.472;
P_amb = 101325;
g = 9.81;

f = @(Pe) ((gamma+1)/2)^(1/(gamma-1))*(Pe/Pcc).^(1/gamma).*((gamma+1)/(gamma-1)*[1-(Pe/Pcc).^((gamma-1)/gamma)]).^(1/2)-1/e_ratio;
%[xvect,it] = bisez(0,100000,1*10^(-9),f);
%Pe = xvect(end);
%options = optimoptions('fzero','FunctionTolerance',1e-9);
[Pe] = fzero(f,50000);%,options); 
Ve = ((2*gamma/(gamma-1))*R/mmol*Tcc*[1-(Pe/Pcc)^((gamma-1)/gamma)])^(1/2);
rho_cc = Pcc/R/Tcc*mmol;
rho_e = rho_cc*(Pe/Pcc)^(1/gamma);
mdot_p = mprop/tburn/nengines;
A_e = mdot_p/rho_e/Ve;
Te = Pe/R/rho_e*mmol;
a_e = (gamma*R/mmol*Te)^(1/2);
Me = Ve/a_e;
IspSL = Ve/g + (Pe - P_amb)/(rho_e*Ve*g);
IspVac = Ve/g + Pe/(rho_e*Ve*g);

out = struct('Pe',      Pe,    ...    
             'Ve',      Ve,    ...    
             'rho_cc',  rho_cc,...    
             'rho_e',   rho_e, ...    
             'mdot_p',  mdot_p,...    
             'A_e',     A_e,   ...    
             'a_e',     a_e,   ...    
             'Te',      Te,    ...    
             'Me',      Me,    ...    
             'IspSL',   IspSL, ...    
             'IspVac',  IspVac);      
