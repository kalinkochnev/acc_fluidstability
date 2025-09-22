function load_envs(enable)
    % Enable environmental variables
    % If this is enabled, if any of the variables are UNSELECTED
    % or missing, then the program will exit.
    %
    % This function is useful when you want to submit jobs to 
    % SLURM but do not want to constantly modify source files
    % before submitting.
    fprintf("(2) Checking for environmental variables...");
    if (enable)
        fprintf("Env file was loaded\n");
        loadenv("../.env")
    else
        fprintf("Env file is selected to not load\n");
        return
    end

    global IS_HPC;
    % if env variables are not available then
    if ~(isenv("IS_HPC"))
        fprintf("IS_HPC is not defined in .env file. Exiting.\n");
        exit;
    end
    IS_HPC = strcmp(getenv("IS_HPC"), "true");

    fprintf("Success!\n");
    fprintf("IS_HPC: %s\n", string(IS_HPC));

    if IS_HPC 
        fprintf("(1b) Adding YALMIP to path...\n");
        if ~(isenv("HOMEDIR"))
            fprintf("HOMEDIR is not specified in .env file\n. Exiting.");
            exit;
        end
        home_dir = getenv("HOMEDIR");

        %add path of YALMIP
        yalmip_path = strcat(home_dir, 'YALMIP-master');
        addpath(genpath(yalmip_path));

        fprintf("Success!\n");
        fprintf("YALMIP path: %s\n", yalmip_path);

        %add path of mosek
        fprintf("(4) Adding MOSEK to path...");
        mosek_path = strcat(home_dir, 'mosek/');
        system(strcat('export PATH=$PATH:', mosek_path, '10.1/tools/platform/linux64x86/bin'));
        addpath(genpath(mosek_path));
        fprintf("Success!\n");

        fprintf("MOSEK path: %s\n", mosek_path);
    end

    global FLUID_SYSTEM;
    global FLUID_SOLVER;

    % Load FLUID_SYSTEM configuration
    disp("(2b) Checking environmental variables are set correctly");
    if (isenv("FLUID_SYSTEM"))
        switch upper(getenv("FLUID_SYSTEM"))
            case 'RADKO_GRL'
                FLUID_SYSTEM = SystemSelect.RADKO_GRL;
                disp("FLUID_SYSTEM=RADKO_GRL");
            case 'RADKO_SHEAR_INSTABILITY'
                FLUID_SYSTEM = SystemSelect.RADKO_SHEAR_INSTABILITY;
                disp("FLUID_SYSTEM=RADKO_SHEAR_INSTABILITY");

            case 'UNSELECTED'
                disp("ERROR: FLUID_SYSTEM is required but unselected.");
                exit;
        end
    else
        disp("ERROR: FLUID_SYSTEM is required but missing from the file.");
        exit;
    end

    % Load FLUID_SOLVER configuration
    if (isenv("FLUID_SOLVER"))
        switch upper(getenv("FLUID_SOLVER"))
            case 'ODE45'
                disp("FLUID_SOLVER=ODE45");
                FLUID_SOLVER = SolverSelect.ODE45;
            case 'CONVEX'
                disp("FLUID_SOLVER=CONVEX");
                FLUID_SOLVER = SolverSelect.CONVEX;
			case 'FLOQUET'
				disp("FLUID_SOLVER=FLOQUET");
				FLUID_SOLVER = SolverSelect.FLOQUET;
            case 'UNSELECTED'
                disp("ERROR: FLUID_SOLVER is required but unselected.");
                exit;

        end
    else
        disp("ERROR: FLUID_SOLVER is required but missing from the file.");
        exit;
    end
end
