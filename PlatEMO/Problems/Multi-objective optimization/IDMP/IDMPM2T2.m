classdef IDMPM2T2 < PROBLEM
% <2019> <multi> <real> <multimodal>
% Imbalanced distance minimization problem with 2 objectives and 2 segments
% a --- 0.4 --- alpha

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, G. G. Yen, Y. Nojima and N. Masuyama, Handling
% Imbalance Between Convergence and Diversity in the Decision Space in
% Evolutionary Multimodal Multi-Objective Optimization, IEEE Transactions
% on Evolutionary Computation, 2019, Early Access,
% DOI: 10.1109/TEVC.2019.2938557.
%------------------------------- Copyright --------------------------------
% Copyright 2018-2019 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    methods
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 2;
            obj.lower    = [-1,-1];
            obj.upper    = [1,1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            a    = 0.4;
            NP   = 2;
            psize = 0.10;
            center = [-0.50,0.50];
            Points = [center - psize; center + psize];
            [N,~]  = size(X);
            PopObj = NaN(N,obj.M);
            for i = 1 : obj.M
                temp = abs(repmat(X(:,1),1,NP) - repmat(Points(i,:),N,1));
                temp(:,1) = temp(:,1) + 100 * (abs(X(:,2)+0.5)).^2;
                temp(:,2) = temp(:,2) + 100 * (abs(X(:,2)-0.5)).^(2-a);
                PopObj(:,i) = min(temp,[],2);
            end
        end
        function R = GetPF(obj)
            R(:,1) = linspace(0,0.2,100)';
            R(:,2) = 0.2 - R(:,1);
        end
    end
end