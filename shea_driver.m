% SHEA_DRIVER driver for shea_compute.
%
%   Useful for parameter space exploration.
%
%   August, 2020, Dani Schmid

% Initialize
clear variables;

% Constants
[yr, myr, km, R, C2K] = shea_constants();

% Config - parameter space variables
Lithology           = {'anorthite_wet', 'anorthite_dry'};
Vel                 = [1:6]/100/yr;
T_sz                = [600, 650];

% Config parts that are overwritten by the above parameter space variables
% config.disp_vel     = 3/100/yr;         % displacement velocity
% config.lithology    = 'anorthite_wet';  % lithology
% config.t_sz         = 600;              % initial temperature in the middle of the shear zone, used to calculate gradient

% Config - user specified model parameters
config.disp_time    = 5*myr;            % duration

config.h_sz         = .5*km;            % shear zone width
config.h_top        = 55*km;         	% model extent above shear zone (1.5 GP)
config.h_bot        = config.h_top;     % model extent below shear zone

config.np           = 2000;          	% numerical resolution - entire model, i.e. shear zone much lower resolution
config.mech_conv    = 1e-4;            	% mechanical convergence criterion
config.plot_freq    = 0;                % plot frequency, 0 = no plot, inf = plot every tstep, myr/10 = example where we plot every 100'000 years
config.plot_agif	= 'test.gif';     	% if we specify a file name here then an animated gif is written

% Loop through parameter space
model_counter	= 0;
for model_litho = 1:length(Lithology)
    % Current lithology
    config.lithology = Lithology{model_litho};
    
    for model_t_sz = 1:length(T_sz)
        % Current initial temperature center shear zone
        config.t_sz = T_sz(model_t_sz);
        
        % Initilize results
        model_counter   = model_counter + 1;
        results(model_counter).lithology    = config.lithology;
        results(model_counter).t_sz         = config.t_sz;
        
        for model_vel = 1:length(Vel)
            % Current displacement velcoity
            config.disp_vel = Vel(model_vel);
            
            disp(['Computing: ', config.lithology, ' @', num2str(config.t_sz), '[C] & ', num2str(config.disp_vel*100*yr), '[cm/yr]']);
            
            t_record = shea_compute(config);
            
            % Record results
            results(model_counter).Vel(model_vel)   = config.disp_vel;
            results(model_counter).T(model_vel)     = max(t_record.T);
            results(model_counter).Tau(model_vel)	= max(abs(t_record.Tau));
        end
    end
end

% Plot Results
LineStyle   = {'--', '--', '-', '-'};
Marker      = {'none', 'd', 'none', 'd'};

h_fig   = figure;
h_ax    = axes(h_fig);
hold(h_ax, 'on');
for i=1:length(results)
    h_p = plot(h_ax, results(i).Vel*100*yr, results(i).T, 'DisplayName', [results(i).lithology, ' @ ', num2str(results(i).t_sz), '[C]']);
    h_p.Color       = [0 0 0];
    h_p.Marker      = Marker{i};
    h_p.LineStyle   = LineStyle{i};
    h_p.LineWidth   = 1;
end
xlabel(h_ax, 'Shear Velocity [cm/yr]');
ylabel(h_ax, 'Max. Temperature Center Shear Zone [C]');
grid(h_ax, 'on');
box(h_ax,  'on');
legend(h_ax, 'Location', 'northeastoutside', 'Interpreter', 'none');

h_fig   = figure;
h_ax    = axes(h_fig);
hold(h_ax, 'on');
for i=1:length(results)
    h_p = plot(h_ax, results(i).Vel*100*yr, results(i).Tau/1e6, 'DisplayName', [results(i).lithology, ' @ ', num2str(results(i).t_sz), '[C]']);
    h_p.Color       = [0 0 0];
    h_p.Marker      = Marker{i};
    h_p.LineStyle   = LineStyle{i};
    h_p.LineWidth   = 1;
end
xlabel(h_ax, 'Shear Velocity [cm/yr]');
ylabel(h_ax, 'Max. Shear Stress [MPa]');
grid(h_ax, 'on');
box(h_ax,  'on');
legend(h_ax, 'Location', 'northeastoutside', 'Interpreter', 'none');
h_ax.YScale = 'log';
