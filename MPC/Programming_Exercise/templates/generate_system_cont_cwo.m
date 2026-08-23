%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ac, Bc] = generate_system_cont_cwo(params)
    % YOUR CODE HERE
    mu = params.model.GravitationalParameter;
    R = params.model.TargetRadius;
    wn = sqrt(mu/R^3);
    Ac = [0 0 0 1 0 0
          0 0 0 0 1 0
          0 0 0 0 0 1
          3* wn^2 0 0 0 2*wn 0
          0 0 0 -2*wn 0 0
          0 0 -wn^2 0 0 0];

    Bc = 1/(params.model.Mass)*[0 0 0
          0 0 0
          0 0 0
          1 0 0
          0 1 0
          0 0 1];
end