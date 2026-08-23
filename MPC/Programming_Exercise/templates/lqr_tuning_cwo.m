%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tuning_struct, i_opt] = lqr_tuning_cwo(x0,Q,params)
    M = size(Q,2);
    R = eye(params.model.nu);

    template = struct( ...
        'InitialCondition', zeros(params.model.nx, 1), ...
        'Qdiag', zeros(params.model.nx, 1), ...
        'MaxAbsPositionXZ', NaN, ...
        'MaxAbsPositionY', NaN, ...
        'MaxAbsThrust', NaN, ...
        'InputCost', NaN, ...
        'MaxFinalPosDiff', NaN, ...
        'MaxFinalVelDiff', NaN, ...
        'TrajFeasible', false);
    tuning_struct = repmat(template, M, 1);

    for i = 1:M
        ctrl = LQR(diag(Q(:,i)), R, params);
        [Xt,Ut] = simulate(x0, ctrl, params);
        [s_max, y_max, u_max, J_u, df_max, vf_max, traj_feas] = traj_constraints_cwo(Xt,Ut,params);

        tuning_struct(i).InitialCondition = x0;
        tuning_struct(i).Qdiag = Q(:,i);
        tuning_struct(i).MaxAbsPositionXZ = s_max;
        tuning_struct(i).MaxAbsPositionY = y_max;
        tuning_struct(i).MaxAbsThrust = u_max;
        tuning_struct(i).InputCost = J_u;
        tuning_struct(i).MaxFinalPosDiff = df_max;
        tuning_struct(i).MaxFinalVelDiff = vf_max;
        tuning_struct(i).TrajFeasible = traj_feas;
    end

    feasible = [tuning_struct.TrajFeasible];
    feasible_costs = [tuning_struct.InputCost];
    feasible_costs(~feasible) = Inf;

    [optimal_cost, i_opt] = min(feasible_costs);
    if isinf(optimal_cost)
        i_opt = NaN;
    end

end
