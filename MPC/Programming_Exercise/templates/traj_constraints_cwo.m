%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s_max, y_max, u_max, J_u, df_max, vf_max, traj_feas] = traj_constraints_cwo(X_t,U_t,params)
    % YOUR CODE HERE
    s_max = max(abs(X_t([1,3], :)), [], 'all');
    y_max = max(abs(X_t(2,:)));
    u_max = max(abs(U_t), [], 'all');

    J_u = sum(U_t.^2, 'all');

    df_max = norm(X_t(1:3, end), 2);
    vf_max = norm(X_t(4:6, end), 2);

    traj_feas = ...
    s_max <= params.constraints.MaxAbsPositionXZ && ...
    y_max <= params.constraints.MaxAbsPositionY && ...
    u_max <= params.constraints.MaxAbsThrust && ...
    df_max <= params.constraints.MaxFinalPosDiff && ...
    vf_max <= params.constraints.MaxFinalVelDiff;

end

