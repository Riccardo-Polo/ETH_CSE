clear;
clc;
%close all;

projectDir = 'C:\Users\ricca\Documents\GitHub\ETH_CSE\MPC\Programming_Exercise';
addpath(genpath(projectDir));

params = generate_params_cwo();

x0 = params.exercise.InitialConditionA;
R  = eye(params.model.nu);

q0 = params.exercise.QdiagOptA;

% Every column is a different vector q
Qdiags = [ ...
    q0, ...
    q0 .* [10; 10; 10; 1; 1; 1], ...   % more weight on position
    q0 .* [1; 1; 1; 10; 10; 10], ...    % more weight on velocity
    [71.6500 ; 0.0469 ; 231.8200 ; 0.0003 ;0.0029 ; 0.0029 ]
];

labels = [
    "Initial value"
    "More weight on position"
    "More weight on velocity"
];


    q = Qdiags(:, 4);
    Q = diag(q);

    ctrl = LQR(Q, R, params);

    [Xt, Ut, ctrl_info] = simulate(x0, ctrl, params);

    [figTime, ~, figPos, ~] = ...
        plot_trajectory_cwo(Xt, Ut, ctrl_info, params);

    figure(figTime);

    figure(figPos);