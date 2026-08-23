%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H_u, h_u, H_x, h_x] = generate_constraints_cwo(params)
    % YOUR CODE HERE
    smax = params.constraints.MaxAbsPositionXZ;
    ymax = params.constraints.MaxAbsPositionY;
    umax = params.constraints.MaxAbsThrust;

    H_x = [1 0 0 0 0 0
           -1 0 0 0 0 0
           0 1 0 0 0 0
           0 -1 0 0 0 0
           0 0 1 0 0 0
           0 0 -1 0 0 0];

    h_x = [smax;smax;ymax;ymax;smax;smax];

    H_u = [1 0 0
           -1 0 0
           0 1 0
           0 -1 0
           0 0 1
           0 0 -1];

    h_u = umax*ones(6,1);
end