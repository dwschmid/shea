function [yr, myr, km, R, C2K] = shea_constants()
% shea_constants   constants required by the shea tools
%
%   August, 2020, Dani Schmid

yr          = 365*24*60*60;         % s
myr         = 1e6*yr;               % s
km          = 1000;                 % m
R           = 8.314;                % J/mol/K, universal gas constant
C2K         = 273.15;               % Celsius to Kelvin conversion