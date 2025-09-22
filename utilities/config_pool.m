function n_samples = config_pool(hpc_cores, hpc_samples, local_cores, local_samples)
   % Change performance settings based on environmental variables
   % Specify number of number of cores, number of samples for 
   % hpc and non-hpc environments.
   % Make sure to call configureDeps() beforehand.
   global IS_HPC;

   pool_exists = ~isempty(gcp('nocreate')); % do not re-create thread pool if it already exists
   if IS_HPC
       n_samples=hpc_samples;
       if ~pool_exists
           parpool(hpc_cores);
       end
   else
       n_samples=local_samples;
       if ~pool_exists
           parpool(local_cores);
       end
   end
end
