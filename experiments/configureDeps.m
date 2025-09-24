function configureDeps()
    % Reload our project dependencies after every execution

    fprintf("--------------- Dependency Configuration ---------------\n");
    fprintf("Current path: %s\n", pwd);
    
    fprintf("(1) Adding utilities and solvers to matlab path...");
    % include utilities from directory above this 
    addpath('../utilities');

    % include solvers from directory above this 
    addpath('../solvers');
    fprintf("Success!\n");
end
