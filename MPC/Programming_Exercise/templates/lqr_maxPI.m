%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, h] = lqr_maxPI(Q,R,params)
	% YOUR CODE HERE
    ctrl = LQR(Q, R, params);
    K = ctrl.K;
    
    A  = params.model.A;
    B  = params.model.B;
    
    Hx = params.constraints.StateMatrix;
    hx = params.constraints.StateRHS;
    
    Hu = params.constraints.InputMatrix;
    hu = params.constraints.InputRHS;
    
    Acl = A + B*K;
    
    % Admissible states must satisfy both:
    % Hx*x <= hx
    % Hu*K*x <= hu
    H_adm = [Hx; Hu*K];
    h_adm = [hx; hu];
    
    X_adm = Polyhedron('A', H_adm, 'b', h_adm);
    
    % Autonomous closed-loop system: no B matrix and no free input
    system = LTISystem('A', Acl);
    
    % Maximum positively invariant set contained in X_adm
    X_LQR = system.invariantSet('X', X_adm);
    
    % Remove redundant inequalities
    X_LQR.minHRep();
    
    % Return the H-representation
    H = X_LQR.A;
    h = X_LQR.b;

end
