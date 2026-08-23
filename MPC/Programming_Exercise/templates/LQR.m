%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026, Marcell Bartos, Sabrina Bodmer, Marco Heim, Institute for Dynamic Systems and Control, ETH Zurich.
%
% All rights reserved.
%
% Please see the LICENSE file that has been included as part of this package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef LQR
    properties
        K
    end
    
    methods
        %constructor
        function obj = LQR(Q,R,params)
            % YOUR CODE HERE
            P = idare(params.model.A , params.model.B , Q , R ,[]);
            obj.K = -((params.model.B)'*P*params.model.B + R)\(params.model.B)'*P*params.model.A; 
        end
        
        function [u, ctrl_info] = eval(obj,x)
            % YOUR CODE HERE
            u = obj.K * x;
            ctrl_info = struct('ctrl_feas',true);
        end
    end
end