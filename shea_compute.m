function t_record = shea_compute(config)
% shea_compute shearzone model with thermo-mechanical feedback in 1D.
%
%   This model consists of two coupled components:
%
%   1) Thermal Model
%      The thermal model is initialized to a steady state ambient state
%      based on a temperature gradient assumption. Transient heat diffusion
%      is then calculated including the effect of shear heating generated
%      by the mechanical deformation.
%
%   2) Mechanical Model
%      The mechanical model assumes that we know the width of the shear
%      zone and that it stays constant. We also assume that we know the
%      shear velocity. The effective viscosity depends locally on the
%      temperature and the shear rate intensity. The internal deformation
%      of the shear zone is solved so that shear stress is constant.
%
%   A self-explanatory example of the input config is given below for the
%   case where no input is specified (nargin==0).
%
%   May, 2020, Dani Schmid

% Constants
[yr, myr, km, R, C2K] = shea_constants();

% Nargin - setup standard config in case none provided
if nargin==0
    config.lithology    = 'anorthite_wet';      % lithology
    config.disp_time    = 5*myr;                % duration
    config.disp_vel     = 5/100/yr;             % displacement velocity
    
    config.h_sz         = 1*km;                 % shear zone width
    config.h_top        = 55*km;                % model extent above shear zone (1.5 GP)
    config.h_bot        = config.h_top;         % model extent below shear zone
    
    config.t_sz         = 600;                  % initial temperature in the middle of the shear zone, used to calculate gradient
    
    config.np           = 2000;                 % numerical resolution - entire model, i.e. shear zone much lower resolution
    config.mech_conv    = 1e-4;                 % mechanical convergence criterion
    config.plot_freq    = myr/10;               % plot frequency, 0 = no plot, inf = plot every tstep
    config.plot_agif	= 'shea_example.gif'; 	% if we specify a file name here then an animated gif is written
end

% Material database
materials	= shea_materials();

% Extract flowlaw parameters
A           = materials{config.lithology,'A'};
n           = materials{config.lithology,'n'};
Q           = materials{config.lithology,'Q'};
f_H2O       = materials{config.lithology,'f_H2O'};
r           = materials{config.lithology,'r'};
rho        	= materials{config.lithology,'rho'};
cp        	= materials{config.lithology,'cp'};
lambda      = materials{config.lithology,'lambda'};

% Diffusivity
D           = lambda/rho/cp;

% Mesh
z_max       = config.h_top + config.h_sz + config.h_bot;
Z           = linspace(0, z_max, config.np)';
dz          = Z(2)-Z(1);

% Shear zone indices
Ind_sz      = find(config.h_top<=Z & Z<=(config.h_top+config.h_sz) );
n_sz        = length(Ind_sz);

% Shear zone indices + part of host rock on each side
% Used for plotting
Ind_sz2     = find((config.h_top-config.h_sz)<=Z & Z<=(config.h_top+2*config.h_sz) );

% Shear zone Z
Z_sz        = Z(Ind_sz);
Z_sz_inter	= (Z_sz(1:end-1)+Z_sz(2:end))/2;  % Interval averages

% Background temperature
t_surf      = 0;
t_grad      = config.t_sz/(config.h_top+config.h_sz/2);
T           = t_surf + Z*t_grad;
T_ori       = T;

% Initial, linear velocity profile across shearzone
Vx          = linspace(0, config.disp_vel, n_sz)';

% Ambient Shear Rate
gamma_r     = config.disp_vel/config.h_sz;

% Sqrt second invariant.
er_ii       = gamma_r/2;

% Conversion factor - Gerya p. 77, eqn. 6.10
F           = 1/( 2^((n-1)/n) * 3^((n+1)/(2*n)));

% Effective A
% Note
% F is excluded from this - like Gerya p. 77
A_eff   = A*f_H2O^r;

%% Time loop
dt_stable	= 0.5*dz^2/D; % stable timestep explicit temperature solver
Time        = linspace(0, config.disp_time, ceil(config.disp_time/dt_stable)+1);
dt          = Time(2)-Time(1);

% Time evolution recorders
t_record.Time	= Time;                            % time
t_record.Tau	= NaN(size(Time));                 % shear stress
t_record.T      = NaN(size(Time));                 % temperature
t_record.T(1)	= T_ori(Ind_sz(round(end/2)));     % read out the middle of the shear zone

for tstep=1:length(Time)-1
    % Mechanics - Shear Zone Only
    T_inter	= (T(Ind_sz(1:end-1))+T(Ind_sz(2:end)))/2; % Interval average temperature
    Vx_old  = Vx;
    niter   = 0;
    while niter==0 || max(Vx-Vx_old)/config.disp_vel>config.mech_conv
        % Iteration counter
        niter   = niter+1;
        
        % Shear rate - interval
        dvxdz   = (Vx(2:end)-Vx(1:end-1))/dz;
        
        % Shear rate invariant - interval
        er_ii   = abs(dvxdz)/2;
        
        % Effective viscosity - interval
        Eta_eff	= F * A_eff.^(-1/n) .* er_ii.^(1/n-1) .* exp(Q/n./R./(T_inter+C2K));
        
        % Setup matrix - shear stress must be constant
        Ind     = 2:n_sz-1;
        A       = ...
            sparse(Ind, Ind-1,  Eta_eff(Ind-1)/dz                  , n_sz, n_sz) + ...
            sparse(Ind, Ind  , -Eta_eff(Ind-1)/dz - Eta_eff(Ind)/dz, n_sz, n_sz) + ...
            sparse(Ind, Ind+1,  Eta_eff(Ind)/dz                    , n_sz, n_sz);
        
        % Rhs intialize
        Rhs         = zeros(n_sz, 1);
        
        % Velocity boundary conditions
        % Note
        % Multiply the entries with a representative viscosity so
        % that the matrix is not singular
        A(1,1)      = 1*A(2,1);
        Rhs(1)      = config.disp_vel*A(2,1);
        
        A(end,end)	= 1*A(2,1);
        Rhs(end)	= 0;
        
        % Solve
        Vx_old  = Vx;
        Vx      = A\Rhs;
    end
    
    % Compute Shear Stress
    dvxdz   = (Vx(2:end)-Vx(1:end-1))/dz;
    er_xz   = dvxdz/2;
    Tau     = 2*Eta_eff.*er_xz;
    
    % Shear heating term - split from intervals into nodes
    SH_sz               = Tau.*er_xz;
    SH                  = zeros(size(T));
    SH(Ind_sz(1:end-1)) = SH(Ind_sz(1:end-1)) + .5*SH_sz;
    SH(Ind_sz(2:end  )) = SH(Ind_sz(2:end  )) + .5*SH_sz;
    
    % Temperature - explicit solve
    % We assume constant conductivity throughout
    T_old = T;
    for i = 2:config.np-1
        T(i) = T_old(i) + dt*D*(T_old(i+1)-2*T_old(i)+T_old(i-1))/dz^2 + dt/(rho*cp)*SH(i);
    end
    
    % Record time evolution
    % Shear stress is constant, just take the first one
    t_record.Tau(tstep+1)	= Tau(1);                   % Shear stress is constant per timestep, does not matter where we read it
    t_record.T(tstep+1)  	= T(Ind_sz(round(end/2)));  % Temperature in the center of the shear zone
    
    % Plot
    if config.plot_freq~=0 && (tstep==1 || mod(tstep+1, round(config.plot_freq/dt))==0 || tstep==length(Time)-1 || isinf(config.plot_freq))
        h_f = figure(1);
        clf(h_f);
        h_f.Color = 'w';
        h_tl = tiledlayout(h_f, 1,6);
        title(h_tl, [config.lithology, ' Vel:', num2str(config.disp_vel*yr),' [m/yr] @', num2str(Time(tstep+1)/1000/yr),' [kyr] Displ:', num2str(config.disp_vel*Time(tstep+1)/1000), ' [km]'], 'interpreter', 'none', 'FontWeight', 'bold');
        
        i = 1;
        h_sp(i) = nexttile;
        plot(h_sp(i), T(Ind_sz2), Z(Ind_sz2)/km);
        hold(h_sp(i), 'on');
        plot(h_sp(i), T(Ind_sz), Z(Ind_sz)/km, '.r');
        h_sp(i).YDir = 'rev';
        title(h_sp(i), 'Temperature', 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        i = i + 1;
        h_sp(i) = nexttile;
        Delta_T     = T-T_ori;
        plot(h_sp(i), Delta_T(Ind_sz2), Z(Ind_sz2)/km);
        hold(h_sp(i), 'on');
        plot(h_sp(i), Delta_T(Ind_sz), Z(Ind_sz)/km, '.r');
        h_sp(i).YDir = 'rev';
        title(h_sp(i), 'Delta Temp', 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        i = i + 1;
        h_sp(i) = nexttile;
        plot(h_sp(i), Vx*yr*100, Z_sz/km);
        h_sp(i).YDir = 'rev';
        h_sp(i).YLim = h_sp(1).YLim;
        title(h_sp(i), 'Velocity [cm/yr]', 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        i = i + 1;
        h_sp(i) = nexttile;
        plot(h_sp(i), Eta_eff, Z_sz_inter/km);
        h_sp(i).YDir = 'rev';
        h_sp(i).YLim = h_sp(1).YLim;
        title(h_sp(i), 'Viscosity [Pas]', 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        i = i + 1;
        h_sp(i) = nexttile;
        plot(h_sp(i), -er_xz, Z_sz_inter/km);
        h_sp(i).YDir = 'rev';
        h_sp(i).YLim = h_sp(1).YLim;
        title(h_sp(i), 'Strain Rate [s^-1]', 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        i = i + 1;
        h_sp(i) = nexttile;
        plot(h_sp(i), -Tau/1e6, Z_sz_inter/km);
        h_sp(i).YDir = 'rev';
        h_sp(i).YLim = h_sp(1).YLim;
        tau_max     = round(max(abs(Tau/1e6)));
        h_sp(i).XLim = [0, tau_max*1.2];
        title(h_sp(i), ['Tau: ', num2str(tau_max), ' MPa'], 'FontWeight', 'normal')
        grid(h_sp(i), 'on');
        
        linkaxes(h_sp, 'y');
        
        % Record animated gif
        if ~isempty(config.plot_agif)
            frame = getframe(h_f);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ~exist(config.plot_agif, 'file')
                imwrite(imind, cm, config.plot_agif, 'gif', 'LoopCount', 1, 'DelayTime', 1);
            else
                imwrite(imind, cm, config.plot_agif, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
            end
        end
    end
end
