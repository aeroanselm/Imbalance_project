function [rho] = EDMA(h)

% --- rho = atmosphere_density(h) ---
%
% Exponentially Decaying Model Atmosphere
%
% INPUT:
%
%       h[1,1] = altitude [km]
%
% OUTPUT:
%
%       rho[1,1] = density [kg/m^3]

% Reference heights [km]

h_0 = ...
[ 0 25 30 35 40 45 ...
  50 55 60 65 70 75 ...
  80 85 90 95 100 110 ...
  120 130 140 150 160 180 ...
  200 250 300 350 400 450 ...
  500 600 700 800 900 1000];

% Reference densities [kg/m^3]

rho_0 = ...
[ 1.225 3.899e-2 1.774e-2 8.279e-3 3.972e-3 1.995e-3 ...
  1.057e-3 5.821e-4 3.206e-4 1.718e-4 8.770e-5 4.178e-5 ...
  1.905e-5 8.337e-6 3.396e-6 1.343e-6 5.297e-7 9.661e-8 ...
  2.438e-8 8.484e-9 3.845e-9 2.070e-9 1.224e-9 5.464e-10 ...
  2.789e-10 7.248e-11 2.418e-11 9.158e-12 3.725e-12 1.585e-12 ...
  6.967e-13 1.454e-13 3.614e-14 1.170e-14 5.245e-15 3.019e-15];

% Scale heights [km]

H = ...
[ 8.44 6.49 6.75 7.07 7.47 7.83 ...
  7.95 7.73 7.29 6.81 6.33 6.00 ...
  5.70 5.41 5.38 5.74 6.15 8.06 ...
  11.6 16.1 20.6 24.6 26.3 33.2 ...
  38.5 46.9 52.5 56.4 59.4 62.2 ...
  65.8 79.0 109.0 164.0 225.0 268.0];

% Density computation

if imag(h) ~= 0 % Avoid numerical errors
    h = abs(h);
end

% Initialization

persistent k;

if isempty(k)
  k = 0;
end

% Calculation

if h >= 1000
    k = 36;
elseif h < 0
    k = 1;
else
    for j = 1:35
        if (h >= h_0(j)) && (h < h_0(j+1))
            k = j;
        end
    end
end

rho = rho_0(k)*exp(-(h - h_0(k))/H(k));