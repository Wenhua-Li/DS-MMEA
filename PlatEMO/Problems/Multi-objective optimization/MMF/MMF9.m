classdef MMF9 < PROBLEM
% <2018> <multi> <real> <multimodal>
% Multi-modal multi-objective test function from Liu/Yen/Gong 2018
% Two PS segments defined by sin(2*pi*x2)^6 = 1 at x2 = 0.25 and 0.75

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong. A Multi-Modal Multi-Objective
% Evolutionary Algorithm Using Two-Archive and Recombination Strategies.
% IEEE Transactions on Evolutionary Computation, 2018, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    methods
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 2;
            obj.lower    = [0.1 0.1];
            obj.upper    = [1.1 1.1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            g  = 2 - sin(2*pi*X(:,2)).^6;
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g ./ X(:,1);
        end
        function R = GetOptimum(obj,N)
            % PF: f2 = 1/f1 on the global PS (g=1)
            f1 = linspace(0.1,1.1,N)';
            R  = [f1, 1./f1];
        end
        function R = GetPF(obj)
            R(:,1) = linspace(0.1,1.1,100)';
            R(:,2) = 1 ./ R(:,1);
        end
    end
end