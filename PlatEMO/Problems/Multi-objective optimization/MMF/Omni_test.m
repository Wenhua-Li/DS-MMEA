classdef Omni_test < PROBLEM
% <2018> <multi> <real> <multimodal>
% Omni-test benchmark (Deb & Tiwari 2005), arbitrary number of variables

%------------------------------- Reference --------------------------------
% K. Deb and S. Tiwari. Omni-optimizer: a procedure for single and
% multi-objective optimization. In Evolutionary Multi-Criterion
% Optimization, Springer, 2005, pp. 47-61.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    methods
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 3;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = 6*ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            [N,n] = size(X);
            PopObj = zeros(N,2);
            for j = 1 : N
                x = X(j,:);
                f = zeros(2,1);
                for i = 1 : n
                    f(1) = f(1) + sin(pi*x(i));
                    f(2) = f(2) + cos(pi*x(i));
                end
                PopObj(j,:) = f';
            end
        end
        function R = GetOptimum(obj,N)
            % Omni-test PF: when all x_i = 0, f1 = 0, f2 = n
            % When all x_i = 0.5 (or any integer), f1 = 0, f2 = n
            R = [0, obj.D];
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end