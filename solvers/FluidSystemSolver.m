% inheriting by handle makes it so it is accessed by reference
classdef FluidSystemSolver < handle
    properties
        params % fluid system parameters to use for the solver
        sysSize % the size of the system matrix to solve for
    end

    methods
        function obj = FluidSystemSolver(fsystem_data) 
            obj.params = fsystem_data;

            switch obj.params.SYSTEM
                case SystemSelect.RADKO_GRL
                    obj.sysSize = 5;
                case SystemSelect.RADKO_SHEAR_INSTABLITY
                    obj.sysSize = 4;
            end
        end

        function M = evaluateSys(obj, M, t, t_ind)

            % updates the matrix M to a certain time step
            param = obj.params;

            m = param.m0 - param.Bu(t)*param.k - param.Bv(t)*param.l; % time dependent vertical wave number
            
            c = param.k^2 + param.l^2 + m^2;
           
            % Choose matrix based on FLUID_SYSTEM configuration
            switch param.SYSTEM
                case SystemSelect.RADKO_GRL
                    M(1, 1, t_ind) = -c;
                    M(1, 2, t_ind) = 0;
                    M(1, 3, t_ind) = 0;
                    M(1, 4, t_ind) = 0;
                    M(1, 5, t_ind) = 1;

                    M(2, 1, t_ind) = 0;
                    M(2, 2, t_ind) = -param.tau * c;
                    M(2, 3, t_ind) = 0;
                    M(2, 4, t_ind) = 0;
                    M(2, 5, t_ind) = param.Rp;

                    M(3, 1, t_ind) = -param.Pr * param.k * m / c;

                    M(3, 2, t_ind) = param.Pr * param.k * m / c;
                    M(3, 3, t_ind) = (param.f * param.k * param.l / c) - param.Pr * c;
                    % for simple case velocity is 0: 
                    M(3, 4, t_ind) = param.f - (param.f * param.k^2 / c);
                    M(3, 5, t_ind) = (param.k^2 * param.Au(t) + param.k * param.l * param.Av(t))/c - param.Au(t);

                    M(4, 1, t_ind) = -param.Pr*param.l*m / c;
                    M(4, 2, t_ind) = param.Pr*param.l*m / c;
                    M(4, 3, t_ind) = (param.f*param.l^2)/c - param.f;
                    M(4, 4, t_ind) = -param.Pr*c - (param.f*param.k*param.l)/c;
                    M(4, 5, t_ind) = (param.k*param.l*param.Au(t) + param.l^2 *param.Av(t))/c - param.Av(t);

                    M(5, 1, t_ind) = -param.Pr*(m^2/c - 1);
                    M(5, 2, t_ind) = param.Pr*(m^2/c - 1);
                    M(5, 3, t_ind) = param.f*param.l*m/c;
                    M(5, 4, t_ind) = -param.f*param.k*m/c;
                    M(5, 5, t_ind) = (param.k*m*param.Au(t) + param.l*m*param.Av(t))/c - param.Pr*c;

                    % in 2D model v=0 so make 4th column zeros
                    M(:, 4) = 0;

                case SystemSelect.RADKO_SHEAR_INSTABLITY
                    M(1, 1, t_ind) = -c;
                    M(1, 2, t_ind) = 0;
                    M(1, 3, t_ind) = 0;
                    M(1, 4, t_ind) = 1;
                    M(2, 1, t_ind) = param.Pr * param.k * m / c;
                    M(2, 2, t_ind) = (param.f * param.k * param.l / c) - param.Pr * c;
                    % for simple case velocity is 0: 
                    M(2, 3, t_ind) = param.f - (param.f * param.k^2 / c);
                    M(2, 4, t_ind) = (param.k^2 * param.Au(t) + param.k * param.l * param.Av(t))/c - param.Au(t);
                    M(3, 1, t_ind) = param.Pr*param.l*m / c;
                    M(3, 2, t_ind) = (param.f*param.l^2)/c - param.f;
                    % for simple case velocity is 0: 
                    M(3, 3, t_ind) = -param.Pr*c + (-param.f*param.k*param.l)/c; %second term should add minus sign
                    M(3, 4, t_ind) = (param.k*param.l*param.Au(t) + param.l^2 *param.Av(t))/c - param.Av(t);
                    M(4, 1, t_ind) = param.Pr*(m^2/c - 1);
                    M(4, 2, t_ind) = param.f*param.l*m/c;
                    % for simple case velocity is 0: 
                    M(4, 3, t_ind) = -param.f*param.k*m/c;
                    M(4, 4, t_ind) = (param.k*m*param.Au(t) + param.l*m*param.Av(t))/c - param.Pr*c;
            end
        end


        function [lambda_opt, P_opt, diagnostics] = optimizeGrowthRate(obj)

            param = obj.params; 
            M = zeros(obj.sysSize, obj.sysSize, param.n_steps);

            t_list=linspace(0,param.T,param.n_steps);
            dt_list = diff(t_list);
            dt = dt_list(1);

            for t_ind=1:length(t_list)
                t = t_list(t_ind);

                % we have to assign to M since matlab does not modify elements in place unless explicitly returned
                M = obj.evaluateSys(M, t, t_ind);
            end

            if norm(M(:,:,end)-M(:,:,1))<1e-10
                M(:,:,end)=M(:,:,1); %if periodic, enforcing exact periodicity
            end

            %%----pre-solving to determine whether system is stable or
            %%unstable by substitute lambda=0 and solving a feasibility
            %%problem. 

            Constraints_lambda0=[];
            for t_ind=1:length(t_list)-1
                P{t_ind} = sdpvar(obj.sysSize, obj.sysSize);
                Constraints_lambda0 = [Constraints_lambda0,P{t_ind}>=eye(obj.sysSize, obj.sysSize)];
            end 
            P{length(t_list)} = P{1};

            for t_ind=1:length(t_list)-1
                %copy the constraint but with lambda=0. If this can be
                %satisifed, then it will be stable and skip the
                %optimization of solving lambda. 
                finite_diff = (P{t_ind + 1} -P{t_ind}) / dt;
                Constraints_lambda0 = [Constraints_lambda0, M(:,:,t_ind)'*P{t_ind}+P{t_ind}* M(:,:,t_ind) + finite_diff <=0]; 
            end
            opts=sdpsettings('solver','mosek','verbose',0);

            diagnostics=optimize(Constraints_lambda0,[],opts); 
            
            if diagnostics.problem==1 %lambda=0 is infeasible, then optimize to solve a positive lambda

                lambda=sdpvar(1,1); %growth rate that is going to be optimized
                
                Objective=lambda; %set the optimization objective. In default, it will minimize this objective. Here it is just lambda
                Constraints=[];
                for t_ind=1:length(t_list)-1
                    %add the constraints for such that Lyapunov inequality is satisfied for each time t. 
    
                    % P(t) matrix that can form a Lyapunov function as V=x^T P(t) x
                    P{t_ind} = sdpvar(obj.sysSize, obj.sysSize);
    
                    Constraints = [Constraints,P{t_ind}>=eye(obj.sysSize, obj.sysSize)]; %the constraints that P is positive definite. Here identity is because the problem is homogeneous (the problem does not change after certain rescaling) so we can rescale P-\epsilon I >=0 with arbitrary \epsilon. So just pick \epsilon=1
                end
                
                P{length(t_list)} = P{1}; % ensure that the P matrix is periodic
                
                for t_ind=1:length(t_list)-1
                    finite_diff = (P{t_ind + 1} -P{t_ind}) / dt;
                    Constraints = [Constraints, M(:,:,t_ind)'*P{t_ind}+P{t_ind}* M(:,:,t_ind) + finite_diff - 2*lambda*P{t_ind}<=0];
                end
                %solve the optimization problem. Here, it is actually a bilinear problem (due to lambda*P term)
                %instead of a standarf semi-definite programming, so bisection is required
                %to optimize over lambda. 
                opts=sdpsettings('solver','mosek');
                opts.bisection.absgaptol=1e-5; %more accuracy in bisection
                diagnostics = bisection(Constraints,Objective, opts);
                   
                %the optimal growth rate that provide an upper bound on the growth rate of
                %solutions; i.e., ||x(t)||<= C e^{\lambda t} ||x(0)||
                lambda_opt=value(lambda);
                
                % return the P matrices for each time step
                for t_ind=1:length(t_list) % we do +1 b/c we added another P to ensure periodic
                    P_opt{t_ind} = value(P{t_ind});
                end
                
                %P_opt=value(P); %the P associated with optimal lambda, may be useful in the future. 
                %------------------------
            else
                lambda_opt=0;
                for t_ind=1:length(t_list) % we do +1 b/c we added another P to ensure periodic
                    P_opt{t_ind} = value(P{t_ind});
                end
            end 
        end

		function [growth_rate] = floquetGrowthRate(obj) 
			param = obj.params;
			t_list=linspace(0,param.T,param.n_steps+1);
			dt=param.T/param.n_steps;

			%compute fundamental solution matrix
			Phi=eye(obj.sysSize, obj.sysSize);
			for t_ind=1:length(t_list)
				t=t_list(t_ind);

				% Have temporary matrix to store
				M = zeros(obj.sysSize, obj.sysSize, 1);
				M = obj.evaluateSys(M, t,1);
				Phi=expm(M*dt)*Phi;
			end
			floq=eig(Phi);%floquet multiplier
			growth_rate=log(max(abs(floq)))/param.T; %energy growth rate
		end


        function [largest_coefficients, largest_t, largest_y, largest_quad_norm, largest_init] = ode45GrowthRate(obj, n_trials)
            % iterate through a bunch of random initial conditions and find the largest growth rate

            M = zeros(obj.sysSize, obj.sysSize, 1);

            function y = F(solver, t, y)
                M = solver.evaluateSys(M, t, 1);
                y = M * y;
            end
            param = obj.params;
            
            % information to keep track of about the largest growth rate
            largest_lambda = -Inf;
            largest_quad_norm = [];
            largest_t = [];
            largest_y = []; % the integration values
            largest_init = zeros(obj.sysSize, 1); % initial conditions for largest growth rate
            largest_coefficients = [];

            for trial=1:n_trials
                % Let the variables in order be [rho, u, v, w] (if sysSize==4)
                % Let variables in order be [T, S, u, v, w] (if sysSize==5)
                a = randn(obj.sysSize, 1);
                b = randn(obj.sysSize, 1);

                init_conds = complex(a, b);
                interval = linspace(0, param.T, param.n_steps); 

                [t,y] = ode45(@(t, y) F(obj, t, y), interval, init_conds);

                assert(obj.sysSize == 4 || obj.sysSize == 5, "The quadratic norm is different depending on the system size. Until a better way to account for these differences has been found this will throw an error");

                quadratic_norm = zeros(length(y), 1);

                for ind=1:length(y)
                    switch obj.params.SYSTEM
                        case SystemSelect.RADKO_GRL
                            u = y(ind, 3);
                            v = y(ind, 4);
                            w = y(ind, 5);
                            T = y(ind, 1);
                            S = y(ind, 2);

                            quadratic_norm(ind) = (abs(u)^2 + abs(v)^2 + abs(w)^2)/4 + param.Pr * abs(S-T)^2 / (4 * (obj.params.Rp - 1));

                        case SystemSelect.RADKO_SHEAR_INSTABLITY
                            u = y(ind, 2);
                            v = y(ind, 3);
                            w = y(ind, 4);
                            rho = y(ind, 1);

                            % (|u|^2 + |v|^2 + |w|^2)/4  + Pr*|rho|^2/4
                            quadratic_norm(ind) = (abs(u)^2 + abs(v)^2 + abs(w)^2)/4 + param.Pr * abs(rho)^2 / 4;
                    end
                end

                % Get coefficients of a line fit through the data.
                % coefficients are in form [m, b]
                coefficients = polyfit(t(obj.params.n_steps/2:end), 0.5* log(quadratic_norm(obj.params.n_steps/2:end)), 1);
                
                lambda = coefficients(1);

                if lambda >= largest_lambda
                    largest_lambda = lambda;
                    largest_init = init_conds;
                    largest_t = t;
                    largest_y = y;
                    largest_quad_norm = quadratic_norm;
                    largest_coefficients = coefficients;
                end

            end
            % largest_lambda
            % largest_init
        end
    
    end
    
end

