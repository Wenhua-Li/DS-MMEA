classdef MMF15_a < PROBLEM
% <2018> <multi> <real> <multimodal>
% Scalable multi-modal multi-objective test function (MMF15-a variant)
% Uses sin-coupling on x_{end-1} and Gaussian-weighted sin^2

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
            t  = -0.5*sin(pi*X(:,end-1)) + X(:,end);
            g  = 2 - exp(-2*log10(2).*((t+1/(2*nP)-0.1)/0.8).^2) .* sin(nP*pi*(t+1/(2*nP))).^2;
            PopObj = repmat(1+g,1,M) .* fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)) ...
                .* [ones(size(g,1),1), sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
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