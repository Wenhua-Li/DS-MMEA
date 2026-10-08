classdef MMF1_z < PROBLEM
% <2018> <multi> <real> <multimodal>
% Multi-modal multi-objective test function (MMF1-z variant)
% Two PS segments with different frequency components

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong. A Multi-Modal Multi-Objective
% Evolutionary Algorithm Using Two-Archive and Recombination Strategies.
% IEEE Transactions on Evolutionary Computation, 2018, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    properties
        POS;
    end
    methods
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 2;
            obj.lower    = [1 -1];
            obj.upper    = [3 1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = abs(X(:,1)-2);
            i1 = X(:,1) <  2;
            i2 = X(:,1) >= 2;
            PopObj(i1,2) = 1 - sqrt(PopObj(i1,1)) + 2*(X(i1,2) - sin(6*pi*PopObj(i1,1)+pi)).^2;
            PopObj(i2,2) = 1 - sqrt(PopObj(i2,1)) + 2*(X(i2,2) - sin(2*pi*PopObj(i2,1)+pi)).^2;
        end
        function R = GetOptimum(obj,N)
            X1 = linspace(1,3,N)';
            obj.POS = [X1, sin(6*pi*abs(X1-2)+pi)];
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - sqrt(R(:,1));
        end
        function R = GetPF(obj)
            R(:,1) = linspace(0,1,100)';
            R(:,2) = 1 - sqrt(R(:,1));
        end
    end
end