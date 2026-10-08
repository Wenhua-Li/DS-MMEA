classdef SYM_PART_rotated < PROBLEM
% <2018> <multi> <real> <multimodal>
% Symmetric part-part benchmark rotated by 45 degrees

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
            obj.lower    = [-20 -20];
            obj.upper    = [20 20];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            a = 1; b = 10; c = 8; w = pi/4;
            [N,~] = size(X);
            PopObj = zeros(N,2);
            for i = 1 : N
                xx1 = cos(w)*X(i,1) - sin(w)*X(i,2);
                xx2 = sin(w)*X(i,1) + cos(w)*X(i,2);
                t1  = sign(xx1) * ceil((abs(xx1) - (a+c/2))/(2*a+c));
                t2  = sign(xx2) * ceil((abs(xx2) - b/2)/b);
                t1  = sign(t1) * min(abs(t1),1);
                t2  = sign(t2) * min(abs(t2),1);
                x1  = xx1 - t1*(c+2*a);
                x2  = xx2 - t2*b;
                PopObj(i,:) = SYM_PART_fun([x1,x2],a);
            end
        end
        function R = GetOptimum(obj,N)
            f1 = 1;
            f2 = linspace(0,3,N)';
            R = [f1*ones(N,1), f2];
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(200);
        end
    end
end

function y = SYM_PART_fun(x,a)
y = zeros(1,2);
y(1) = (x(1)+a)^2 + x(2)^2;
y(2) = (x(1)-a)^2 + x(2)^2;
end