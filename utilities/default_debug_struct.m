function debug_struct = default_debug_struct(fparams)
   % We need to define the layout of the struct so matlab can save
   % debug struct arrays properly.
   % This function provides a default.
   switch fparams.SOLVER
       case SolverSelect.ODE45
           debug_struct = struct('initial', [], 't', [], 'y', [], 'quadratic_norms', [], 'coefficients', []);
       case SolverSelect.CONVEX
           debug_struct = struct('lambda', NaN, 'Popt', [], 'diagnostics', []);

       case SolverSelect.FLOQUET
           debug_struct = struct('lambda', NaN);
   end

end
