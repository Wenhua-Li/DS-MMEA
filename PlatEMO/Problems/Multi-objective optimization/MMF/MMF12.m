classdef MMF12 < PROBLEM
% <2018> <multi> <real> <multimodal>
% Multi-modal multi-objective test function from Liu/Yen/Gong 2018
% PS: 2 segments + nonlinear PF with modulation h(f1, g, q)

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
            obj.lower    = [0 0];
            obj.upper    = [1 1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            q      = 4;
            Alfa   = 2;
            g      = 2 - (sin(2*pi*X(:,2)).^6).*(exp(-2*log10(2).*((X(:,2)-0.1)/0.8).^2));
            PopObj(:,1) = X(:,1);
            h      = 1 - (PopObj(:,1)./g).^Alfa - (PopObj(:,1)./g).*sin(2*pi*q*PopObj(:,1));
            PopObj(:,2) = g.*h;
        end
        function R = GetOptimum(obj,N)
            % When g = 1 (at x2 = 0.25, 0.75), the PF is the curve f2 = 1 - f1^2 - f1*sin(2*pi*q*f1)
            % sampled at g=1.
            f1 = linspace(0,1,N)';
            g  = 1;
            h  = 1 - (f1./g).^2 - (f1./g).*sin(2*pi*4*f1);
            R  = [f1, g.*h];
            % Keep only feasible (PF is upper boundary of h)
            R(R(:,2) < 0,:) = [];
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(200);
        end
    end
end