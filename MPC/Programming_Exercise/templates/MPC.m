%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MPC
    properties
        yalmip_optimizer
    end

    methods
        function obj = MPC(Q,R,N,params)
            nu = params.model.nu;
            nx = params.model.nx;

            A = params.model.A;
            B = params.model.B;
            Hx = params.constraints.StateMatrix;
            hx = params.constraints.StateRHS;
            Hu = params.constraints.InputMatrix;
            hu = params.constraints.InputRHS;

            % define optimization variables
            U = sdpvar(repmat(nu,1,N),ones(1,N),'full');
            X0 = sdpvar(nx,1,'full');

            % Infinite-horizon LQR terminal cost l_f(x) = x'*P*x
            P = idare(A, B, Q, R, []);
            P = (P + P')/2;

            constraints = [];
            objective = 0;
            x = X0;

            for i = 1:N
                % Stage cost for x_{i-1}, u_{i-1}
                objective = objective ...
                    + x' * Q * x ...
                    + U{i}' * R * U{i};

                % State and input constraints at the current stage
                constraints = [constraints, ...
                    Hx * x <= hx, ...
                    Hu * U{i} <= hu];

                % Predicted successor state
                x = A * x + B * U{i};
            end

            % Terminal state constraint and terminal LQR cost
            constraints = [constraints, Hx * x <= hx];
            objective = objective + x' * P * x;

            opts = sdpsettings('verbose',0,'solver','quadprog');
            obj.yalmip_optimizer = optimizer(constraints,objective,opts,X0,{U{1} objective});
        end

        function [u, ctrl_info] = eval(obj,x)
            %% evaluate control action by solving MPC problem, e.g.
            tic;
            [optimizer_out,errorcode] = obj.yalmip_optimizer(x);
            solvetime = toc;
            
            [u, objective] = optimizer_out{:};

            feasible = true;
            if (errorcode ~= 0)
                feasible = false;
            end

            ctrl_info = struct('ctrl_feas',feasible,'objective',objective,'solvetime',solvetime);
        end
    end
end
