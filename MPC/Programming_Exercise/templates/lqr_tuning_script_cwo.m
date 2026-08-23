%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

params = generate_params_cwo();
x0 = params.exercise.InitialConditionA;

% A feasible low-cost point obtained by refining a search around QdiagOptA.
q_base = [71.65; 0.0469; 231.82; 2.8626e-4; 2.8626e-3; 2.8626e-3];

% Search the coupled in-plane weights and the decoupled z weights separately.
position_xy_scale = [1.00, 1.10, 1.25];
velocity_xy_scale = [1.00, 1.10, 1.25];
position_z_scale = [1.00, 1.10, 1.25];
velocity_z_scale = [1.00, 1.10, 1.25];

[position_xy_grid, velocity_xy_grid, position_z_grid, velocity_z_grid] = ...
    ndgrid(position_xy_scale, velocity_xy_scale, position_z_scale, velocity_z_scale);

Q = q_base .* [ ...
    position_xy_grid(:)'; ...
    position_xy_grid(:)'; ...
    position_z_grid(:)'; ...
    velocity_xy_grid(:)'; ...
    velocity_xy_grid(:)'; ...
    velocity_z_grid(:)' ...
];

% Include the reference weights supplied with the project in the study.
Q(:,end+1) = params.exercise.QdiagOptA;

[tuning_struct, i_opt] = lqr_tuning_cwo(x0, Q, params);

assert(~isnan(i_opt), 'No feasible LQR controller was found.');
q = tuning_struct(i_opt).Qdiag;
assert(tuning_struct(i_opt).TrajFeasible, 'The selected trajectory is infeasible.');
assert(tuning_struct(i_opt).InputCost <= 8, ...
    'The selected controller does not satisfy J_u <= 8.');

fprintf('Best candidate: %d of %d\n', i_opt, size(Q,2));
fprintf('Input cost: %.6f\n', tuning_struct(i_opt).InputCost);
fprintf('Maximum thrust: %.6f N\n', tuning_struct(i_opt).MaxAbsThrust);
disp('Selected q:');
disp(q);


%% Save results
current_folder = fileparts(which(mfilename));
save(fullfile(current_folder, "lqr_tuning_script_cwo.mat"), 'q','tuning_struct');
