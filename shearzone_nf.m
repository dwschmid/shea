% SHEARZONE_NF simplified shearzone model with no thermo-mechanical feedback.
%
%   Here the temperature is assumed constant throughout a shear zone of
%   fixed width and the boundary condition is a constant shear velocity.
%
%   Can be used to study eff. viscosity and shear stress as a function of
%   temperature, shear zone width, shear velocity, and rheology.
%
%   August, 2020, Dani Schmid

% Constants
[yr, myr, km, R, C2K] = shea_constants();

% User input
lithology   = 'anorthite_dry';      % 'anorthite_wet' 'anorthite_dry'
h_sz        = 5*km;                 % shear zone width
vel_shear	= 0.01/yr;              % shear velocity [m/s]
T_amb       = linspace(600, 700);   % ambient temperatures in Celsius

% Material database
materials	= shea_materials();

% Extract flowlaw parameters
A           = materials{lithology,'A'};
n           = materials{lithology,'n'};
Q           = materials{lithology,'Q'};
f_H2O       = materials{lithology,'f_H2O'};
r           = materials{lithology,'r'};

% Geometry conversion factor - Gerya (2010), p. 77, eqn. 6.10
F           = 1/( 2^((n-1)/n) * 3^((n+1)/(2*n)));

% Ambient Shear Rate
% Since we don't consider varying temperature this is constant
gamma_r     = vel_shear/h_sz;

% Shear rate invariant
er_ii       = gamma_r/2;

% Effective A
% Note
%  F is excluded from this, i.e. like Gerya (2010), but different to
%  Kiss et al. (2019)
A_eff       = A*f_H2O^r;

% Effective viscosity
eta_eff     = F * A_eff.^(-1/n) .* er_ii.^(1/n-1) .* exp(Q/n./R./(T_amb+C2K));

% Shear stress
tau_xy      = eta_eff*gamma_r;

% Plot
h_f = figure;

h_ax1 = subplot(1,2,1);
plot(h_ax1, T_amb, log10(eta_eff));
xlabel(h_ax1, 'Temperature [Celcius]');
ylabel(h_ax1, 'Viscosity [log10(Pas)])');
title(h_ax1, lithology, 'interpreter', 'none');
grid(h_ax1, 'on');

h_ax2 = subplot(1,2,2);
plot(h_ax2, T_amb, tau_xy);
xlabel(h_ax2, 'Temperature [Celcius]');
ylabel(h_ax2, 'Shear Stress [Pa]');
title(h_ax2, lithology, 'interpreter', 'none');
grid(h_ax2, 'on');
h_ax2.YAxis.Exponent = 6;
