clear all;
close all;
clc;

% Declare global variables

configureDeps();

rng(1001)

% Parameters
omega = 0.5; % Frequency parameter
fparams = FluidSystemData(50*2*pi/omega, 2000);

% Transfer environmental variables to simulation parameters
fparams.IS_HPC = false;
fparams.SYSTEM = SystemSelect.RADKO_GRL;
fparams.SOLVER = SolverSelect.ODE45;

% Set parameters
fparams = experiment1_common(fparams);
fparams.omega = omega;
fparams.m0 = -3.01E-3;
fparams.k = 0.179;
solver = FluidSystemSolver(fparams);

n_random_inits = 50;
[coefficients, t, y_vec, quadratic_norms, init_conds] = solver.ode45GrowthRate(n_random_inits);

% If the growth rate is negative/stable then return log10lambda as -Inf
if real(coefficients(1)) <= 0
	log10lambda = -Inf;
else
	log10lambda = log10(real(coefficients(1)));
end

fprintf("ODE45 (k=%.3f, m0=%.3f): %f\n", fparams.k, fparams.m0, log10lambda);

% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(t), max(t), length(t));

% % Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
grid on;


data{1}.x=t;
data{1}.y=0.5*log(quadratic_norms);
data{2}.x=xFit;
data{2}.y=yFit;
plot_config.xlim_list=[1,0,1000];
plot_config.label_list={1,'$t$','$ln(e)/2$'};
plot_config.linewidth=6;
plot_line(data,plot_config);
