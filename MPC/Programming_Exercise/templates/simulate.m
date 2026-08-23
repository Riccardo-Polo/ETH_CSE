%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xt,Ut,ctrl_info] = simulate(x0, ctrl, params)

% YOUR CODE HERE
% Hint: you can access the control command with ctrl.eval(x(:,i))
Nt = params.exercise.SimHorizon;

Xt = zeros(params.model.nx , params.exercise.SimHorizon + 1 );
Ut = zeros(params.model.nu , params.exercise.SimHorizon);
Xt(:,1) = x0;

for i = 1 : params.exercise.SimHorizon
    [Ut(:,i),info_i] = ctrl.eval(Xt(:,i));

    if i == 1
        ctrl_info = repmat(info_i , 1 , Nt);
    end

    ctrl_info(i) = info_i;
    Xt(:,i+1) =params.model.A * Xt(:,i) + params.model.B * Ut(:,i);
end

end
