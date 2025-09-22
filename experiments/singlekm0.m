% this file is useful when you want to test a single k,m0 pair.
close all;
clear all;
clc;

fprintf("---------- Single optimization condition -----------\n");


% Declare global variables
global IS_HPC; 
global FLUID_SOLVER;
global FLUID_SYSTEM;

% load_env_file = false; % Set to true if using an env file
configureDeps();
% if ~load_env_file
%     IS_HPC = false;
%     FLUID_SOLVER = SolverSelect.CONVEX;
%     FLUID_SYSTEM = SystemSelect.RADKO_GRL;
% end
%load_envs(load_env_file);

rng(1001)

% Parameters
omega = 0.5; % Frequency parameter

% R

n_step = 2000;
fprintf("Testing FluidSystemData with parameter value: %d\n", n_step); 
fparams = FluidSystemData(2*pi/omega, n_step);

% Transfer environmental variables to simulation parameters
fparams.IS_HPC = false;
fparams.SYSTEM = SystemSelect.RADKO_GRL;
fparams.SOLVER = SolverSelect.FLOQUET;

% Set parameters
fparams = experiment1_common(fparams);
fparams.omega = omega;
fparams.m0 = 0;
fparams.k = -0.196;


solver = FluidSystemSolver(fparams);
[lambda_opt, P_opt, diagnostics] = solver.optimizeGrowthRate();
lambda_opt_cvx = log10(real(lambda_opt));


param = fparams; 
M = zeros(solver.sysSize, solver.sysSize, n_step);

t_list=linspace(0,param.T,param.n_steps);
dt_list = diff(t_list);
dt = dt_list(1);

for t_ind=1:length(t_list)
	t = t_list(t_ind);

	% we have to assign to M since matlab does not modify elements in place unless explicitly returned
	M = solver.evaluateSys(M, t, t_ind);
end

results.n_step = n_step;
results.k = fparams.k;
results.m0 = fparams.m0;
results.omega = omega;
results.P_opt = P_opt;
results.lambda_opt = lambda_opt_cvx;
results.At = M;

save("../results/convex.mat", 'results');
fprintf("lambda_opt = %f\n", lambda_opt_cvx);
