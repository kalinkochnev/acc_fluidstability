
function lambda_mat=heatmapGrowthRate(fparams, debug_struct, growthRateFunc, x_vec, y_vec, filename)
    % Iterate through different x and y values and calculate the growth rate.
    % Saves the lambda matrix to <filename>.mat, <filename>.fig, and <filename>.png
    % Debug data is saved as one file per `x_ind` containing debug information for all
    % `y_samples`. File format is "/results/debug/<filename>_debugxind_<xindex>.mat"
    %
    % Parameters:
    %   fparams - initial fluid system parameters
    %   growthRateFunc(fparams, x, y) -> [debug, lambda]: Can modify FluidSystemData 
    %   object (fparams) and run a particular solver. 
    %       `debug` must be a structure with fields consistent across invocations.
    %               it is used to record debugging data for that particular (x, y)
    %       `lambda` is the largest growth rate for this `x` and `y`
    %   x_vec - vector of x-axis values to vary "x" parameter with
    %   y_vec - vector of y-axis values to vary "y" paremeter with


    fprintf("--------------- Begin Heatmap Generation ---------------\n");

    x_samples = length(x_vec);
    y_samples = length(y_vec);

    lambda_mat = zeros(y_samples, x_samples);
    [X,Y]=meshgrid(x_vec, y_vec);

    debug_dir = "../results/debug/";
    if ~exist(debug_dir)
        mkdir(debug_dir);
    end

    parfor x_ind=1:x_samples
        % make a copy of the params to avoid overwriting original values
        copy_fparams = fparams;
        tmpGrowthRateFunc = growthRateFunc;
        tmp_x_vec = x_vec;
        tmp_y_vec = y_vec;

        debug_structs = debug_struct;

        for y_ind=1:y_samples
            fprintf("Progress (%d, %d) ", x_ind, y_ind); 
            % Modifies fparams in order to vary and plot whatever system attributes we want onto the heatmap
            [debug, lambda] = tmpGrowthRateFunc(copy_fparams, tmp_x_vec(x_ind), tmp_y_vec(y_ind));

            % Save the debug information struct into an array of structs
            debug_structs(y_ind) = debug;
            lambda_mat(y_ind, x_ind) = lambda;

            % Periodically save the data in 10% increments
            if mod(y_ind, round(0.1 * x_samples)) == 0
                parsave(strcat(debug_dir, filename, "_debugxind_", string(x_ind), '.mat'), debug_structs);
            end
        end

        % Save matrix with debug data into file of format "/results/debug/<filename>_debugxind_<x index>.mat"
        % In order for matlab to save, need to save into a struct first.
        % See here: https://www.mathworks.com/help/parallel-computing/save-variables-in-parfor-loop.html
        
        parsave(strcat(debug_dir, filename, "_debugxind_", string(x_ind), '.mat'), debug_structs);
    end
    save(strcat('../results/', filename, '.mat'),'lambda_mat');
    surf(X,Y,lambda_mat);
    colormap('jet');
    %zlim([0, 0.02]);

    % Only view if not HPC
    if fparams.IS_HPC
    else
        view(2);
    end
    
    % Save matlab figure file and raw png for convenience
    savefig(strcat('../results/', filename, '.fig'))
    saveas(gcf,strcat('../results/', filename', '.png')) 
end

% This is necessary in order to avoid variable transparency errors with
% matlab partfor. Definitely BS.
% https://stackoverflow.com/a/25293356
function parsave(fname, var_to_save) 
    save(fname, 'var_to_save')
end

