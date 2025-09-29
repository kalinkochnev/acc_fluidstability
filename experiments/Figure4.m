function convergence_experiment(k, m0, n_samples_start, n_samples_end, n_samples_step)

    close all;
    clc;

    % Declare global variables
    global IS_HPC; 
    global FLUID_SOLVER;
    global FLUID_SYSTEM;

    load_env_file = false; % Set to true if using an env file
    configureDeps();
    if ~load_env_file
        IS_HPC = false;
        FLUID_SOLVER = SolverSelect.CONVEX;
        FLUID_SYSTEM = SystemSelect.RADKO_GRL;
    end
    load_envs(load_env_file);

    rng(1001)

    % Parameters
    omega = 0.5; % Frequency parameter

    % Range
    param_range = n_samples_start:n_samples_step:n_samples_end;
    num_params = length(param_range);

    cvx_results = zeros(num_params, 2); 
    floquet_results = zeros(num_params, 2);

    % Start parallel pool
    if isempty(gcp('nocreate'))
        parpool(2); 
    end

    parfor i = 1:num_params
        param_value = param_range(i);
        fprintf("Testing FluidSystemData with parameter value: %d\n", param_value);

        fparams = FluidSystemData(2*pi/omega, param_value);

        % Transfer environmental variables to simulation parameters
        fparams.IS_HPC = false;
        fparams.SYSTEM = SystemSelect.RADKO_GRL;

        % Set parameters
        fparams = experiment1_common(fparams);
        fparams.omega = omega;
        fparams.m0 = m0;
        fparams.k = k;

        % Solve system for convex
        fparams.SOLVER = SolverSelect.CONVEX;
        solver = FluidSystemSolver(fparams);
        [lambda_opt, ~, ~] = solver.optimizeGrowthRate();
        lambda_opt_cvx = real(lambda_opt);
        cvx_results(i, :) = [param_value, lambda_opt_cvx];

        % Solve system for floquet
        fparams.SOLVER = SolverSelect.FLOQUET;
        solver = FluidSystemSolver(fparams);
        lambda_opt = solver.floquetGrowthRate();
        lambda_opt_floquet = lambda_opt;
        floquet_results(i, :) = [param_value, lambda_opt_floquet];

        fprintf("(lambda_opt cvx, lambda_opt floquet) = (%f, %f)\n", lambda_opt_cvx, lambda_opt_floquet);
    end

    % Solve with ode45 for large number of samples
    disp("Evaluating numerical simulation...");
    n_random_inits = 50;
    fparams = FluidSystemData(50*2*pi/omega, 2000); % For ODE45
    fparams = experiment1_common(fparams);
    fparams.omega = omega;
    fparams.m0 = m0;
    fparams.k = k;
    fparams.SOLVER = SolverSelect.ODE45;
    fparams.SYSTEM = SystemSelect.RADKO_GRL;
    solver = FluidSystemSolver(fparams);
    [coefficients, ~, ~, ~, ~] = solver.ode45GrowthRate(n_random_inits);

    if real(coefficients(1)) <= 0
        ode45_lambda = -Inf;
    else
        ode45_lambda = real(coefficients(1));
    end

    % Save results into a struct
    results = struct();
    results.parameters = struct('k', k, 'm0', m0, 'omega', omega, ...
                                'n_samples_start', n_samples_start, ...
                                'n_samples_end', n_samples_end, ...
                                'n_samples_step', n_samples_step, ...
                                'n_random_inits', n_random_inits);
    results.cvx_results = cvx_results;
    results.floquet_results = floquet_results;
    results.ode45_lambda = ode45_lambda;

    save("convergence.mat", 'results');
end

convergence_experiment(0.179, 3.01E-3, 2, 40, 5);
