classdef MMF14 < PROBLEM
% <2018> <multi> <real> <multimodal>
% Scalable multi-modal multi-objective test function
% M and D both set to 3 in the original; can be re-set to higher M/D

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong. A Multi-Modal Multi-Objective
% Evolutionary Algorithm Using Two-Archive and Recombination Strategies.
% IEEE Transactions on Evolutionary Computation, 2018, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    methods
        function Setting(obj)
            obj.M        = 3;
            obj.D        = 3;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            M  = obj.M;
            nP = 2;
            g  = 2 - sin(nP*pi*X(:,end)).^2;
            PopObj = repmat(1+g,1,M) .* fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)) ...
                .* [ones(size(g,1),1), sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
            % Sphere-segment PF on (M-1) sphere of radius 2 at origin offset
            M  = obj.M;
            n1 = ceil(N/M);
            f1 = linspace(0,1,n1)';
            f2 = linspace(0,1,n1)';
            if M == 2
                R = [f1, 1-f1; (1-f1), f1];
            else
                R = [f1, 1-f1, 0.5*ones(n1,1);
                     (1-f1), f1, 0.5*ones(n1,1);
                     0.5*ones(n1,1), 1-f1, f1];
            end
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(300);
        end
    end
end