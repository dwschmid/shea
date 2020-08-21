function props = shea_materials()
% shea_materials  database of visous flow law and thermal parameters
%
%   Only Rybacki & Dresen (2000) anorthite is implemented.
%
%   Note
%   Flow laws describe stress and strain rate relationship in uniaxial
%   compression experiments. In order to convert them into strain
%   rate dependent invariant forms we need to introduce a geometry factor,
%   which is according to Gerya (2010, p. 77, eqn. 6.10)
%
%   F = 1/( 2^((n-1)/n) * 3^((n+1)/(2*n)))
%
%   This expression is different to the one used in Kiss et al. (2019)
%
%   F = 2^(n-1)*3^((n+1)/2).
%
%   This is because in the effective viscosity formulation of Gerya F is
%   kept outside A_eff, while Kiss et al. multiply F into A_eff. See
%   shearzone_nf for an example.
%
%   August, 2020, Dani Schmid

lithology   = {'anorthite_wet'; 'anorthite_dry'};  %
A           = [       3.98e-16;         5.01e-6];  % Pa^(−n-r)*s^−1
n           = [              3;               3];  %
f_H2O       = [              0;               0];  % Pa
r           = [              0;               0];  %
Q           = [         3.56e5;          6.48e5];  % J*mol^-1
lambda      = [            2.2;             2.2];  % W*K^-1*m^-1
rho         = [           2730;            2730];  % kg*m^-3
cp          = [           1050;            1050];  % J*kg^-1

props       = table(A, n, f_H2O, r, Q, lambda, rho, cp, 'RowNames', lithology);