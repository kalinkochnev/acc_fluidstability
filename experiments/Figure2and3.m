% this file can generature figure 2 and 3. Figure 3 can be generated if you
% call heatmapGrowthRate(). Figure 2 can be generated if you call
% ode45Singlekm0Plot().
close all;
clear all;
clc;

fprintf("---------- Thermohaline Figure 1A Experiment -----------\n");

% need to declare usage of global variable
global IS_HPC; 
global FLUID_SOLVER;
global FLUID_SYSTEM;

load_env_file = false; % IMPORTANT! Specify as true if you want to use the env file configuration 
configureDeps();
if ~load_env_file
    IS_HPC = false;
    FLUID_SOLVER = SolverSelect.FLOQUET;
    FLUID_SYSTEM = SystemSelect.RADKO_GRL;
end
load_envs(load_env_file);

fprintf("-------------- Experiment Initialization ---------------\n");

% Order of arguments: (# cores for HPC, # samples in wavenumber for HPC, # cores local development, # samples in wavenumber local)
n_samples = config_pool(126, 126, 2, 126);

% If SLURM_ARRAY_TASK_ID is defined, use that to determine number of sampled time_steps
%n_samples = 1;

% Configure fluid system data and growth rate func based on solver choice
omega = 0.5;
output_file_name = "thermohaline_fig1a_";

switch FLUID_SOLVER
    case SolverSelect.ODE45

        if IS_HPC
            fparams = FluidSystemData(50*2*pi/omega, 2000); % For ODE45
        else
            fparams = FluidSystemData(50*2*pi/omega, 2000); % For ODE45
        end

        growthRateFunc = @ode45GrowthRate;
        output_file_name = strcat(output_file_name, "ode45");

    case SolverSelect.CONVEX

		if isenv("SLURM_ARRAY_TASK_ID")
			time_step_options = [200, 400, 600, 800];
			task_id = getenv("SLURM_ARRAY_TASK_ID");
			fprintf("Task ID#%s Specified!\n", task_id);
			fprintf("Time step options: \n", task_id);
			disp(time_step_options);

			n_time_steps = time_step_options(str2num(task_id));
			fprintf("Chosen time step: %d\n", n_time_steps);
			fparams = FluidSystemData(2*pi/omega, n_time_steps); % For for convex optimization
		else
			fprintf("Task ID not found. Proceeding as normal.\n");
			fparams = FluidSystemData(2*pi/omega, 500); % For for convex optimization
		end
        growthRateFunc = @optimizeGrowthRate;
        output_file_name = strcat(output_file_name, "convex");
	case SolverSelect.FLOQUET
        fparams = FluidSystemData(2*pi/omega, 2000); % For for convex optimization
        growthRateFunc = @floquetGrowthRate;
        output_file_name = strcat(output_file_name, "floquet");
end

% Setup all common values across experiments
fparams = experiment1_common(fparams);

% Transfer environmental variables to simulation parameters
% since global variables are not thread safe.
fparams.IS_HPC = IS_HPC;
fparams.SOLVER = FLUID_SOLVER;
fparams.SYSTEM = FLUID_SYSTEM;

k = linspace(-0.5, 0.5, n_samples);
m0 = linspace(0, 1.5, n_samples); % this is m0
fprintf("Varying parameters: k=(%d,%d) m0=(%d,%d)\n", k(1), k(end), m0(1), m0(end));
fprintf("Size: %d x %d\n", length(k), length(m0));

fparams.jsonencode()
% Plot the heatmap for varying initial conditions
lambda_mat=heatmapGrowthRate(fparams, default_debug_struct(fparams), growthRateFunc, k, m0, output_file_name);
%ode45Singlekm0Plot(fparams, 0.179, 3.01E-3);

% Plot for a single initial condition using ode45
function ode45Singlekm0Plot(params, k, m0)
    [debug, log10lambda] = ode45GrowthRate(params, k, m0)
    debug
    t = debug.t;

    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(t), max(t), length(t));

    % % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(debug.coefficients, xFit);
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
    grid on;


    data{1}.x=debug.t;
    data{1}.y=0.5*log(debug.quadratic_norms);
    data{2}.x=xFit;
    data{2}.y=yFit;
    plot_config.xlim_list=[1,0,1000];
    plot_config.label_list={1,'$t$','$ln(e)/2$'};
    plot_config.linewidth=6;
    plot_line(data,plot_config);
end

function [debug, log10lambda] = ode45GrowthRate(params, x, y) 
    params.k = x;
    params.m0 = y;
    solver = FluidSystemSolver(params);

    % Get a rough estimation when using a personal computer
    if params.IS_HPC
        n_random_inits = 50;
    else
        n_random_inits = 50;
    end

    [coefficients, t, y_vec, quadratic_norms,init_conds] = solver.ode45GrowthRate(n_random_inits);

    % Store any debug information into a struct to be saved
    % The struct field members must remain the same across iterations
    % The actual types/values of the struct can be different.
    debug.initial = init_conds;
    debug.t = t;
    debug.y = y_vec;
    debug.quadratic_norms = quadratic_norms;
    debug.coefficients = coefficients;

    % If the growth rate is negative/stable then return log10lambda as -Inf
    if real(coefficients(1)) <= 0
        log10lambda = -Inf;
    else
        log10lambda = log10(real(coefficients(1)));

    end

	fprintf("ODE45 (k=%.3f, m0=%.3f): %f\n", x, y, log10lambda);

end

function [debug, log10lambda] = optimizeGrowthRate(params, x, y)
    params.k = x;
    params.m0 = y;

    solver = FluidSystemSolver(params);
    [lambda, P_opt, diagnostics] = solver.optimizeGrowthRate();

    debug.lambda = lambda;
    debug.Popt = P_opt;
    debug.diagnostics = diagnostics;

    % discard the imaginary part?
    log10lambda = log10(real(lambda));

    fprintf("Optimize (k=%f, m0=%f): %.4f\n", x, y, log10lambda);
end

function [debug, log10lambda] = floquetGrowthRate(params, x, y) 
    params.k = x;
    params.m0 = y;

    solver = FluidSystemSolver(params);
    [lambda] = solver.floquetGrowthRate();


    debug.lambda = lambda;

    % discard the imaginary part?
    log10lambda = log10(lambda);

    fprintf("Floquet (k=%f, m0=%f): %.4f\n", x, y, log10lambda);
end


