classdef FluidSystemData
    properties
        T % stop time
        n_steps % interval size for time discretization discretization
        Ri % mean richardson number
        Pr % prandtl number
        k % x component wave number
        l % y component wave number
        f
        omega % background shear frequency
        m0 % initial vertical wave number 
        Au 
        Av
        Bu
        Bv
        Rp % background density ratio (from Thermohaline-Shear Instability paper)
        tau % diffusivity ratio (from Thermohaline-Shear Instability paper)

        %%% Experimental 

        %%% ENVIRONMENTAL VARIABLES
        IS_HPC
        SOLVER
        SYSTEM
    end


    methods

        function obj = FluidSystemData(stop_time, n_steps) 
            % stop_time - which time to solve up to
            % n_steps - interval size for discretization
            obj.T = stop_time; % usually 2pi/omega for optimizer, 100*2pi/omega for ode45
            obj.n_steps = n_steps; % by default n_steps = stop_time

            % default Bu Bv is the matlab integrate function
            obj.Bu = @(t) integral(obj.Au, 0, t);
            obj.Bv = @(t) integral(obj.Av, 0, t);
        end

    function json = jsonencode(obj, varargin)
        s=struct("T",obj.T,"n_steps",obj.n_steps,"Ri",obj.Ri,"Pr",obj.Pr,"k",obj.k,"l",obj.l,"f",obj.f,"omega",obj.omega,"m0",obj.m0,"Au",func2str(obj.Au),"Av",func2str(obj.Av),"Bu",func2str(obj.Bu),"Bv",func2str(obj.Bv),"Rp",obj.Rp,"tau",obj.tau,"solver",obj.SOLVER,"system",obj.SYSTEM);

        json = jsonencode(s, "PrettyPrint", true);

        fid = fopen('params.json','w');
        fprintf(fid, '%s', json);
        fclose(fid)
    end
    end
end
