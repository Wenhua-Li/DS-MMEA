classdef IDMPM2T3_e < PROBLEM
% <2019> <multi> <real> <multimodal>
% Imbalanced distance minimization problem with 2 objectives and 2 segments

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, G. G. Yen, Y. Nojima and N. Masuyama, Handling
% Imbalance Between Convergence and Diversity in the Decision Space in
% Evolutionary Multimodal Multi-Objective Optimization, IEEE Transactions
% on Evolutionary Computation, 2019, Early Access,
% DOI: 10.1109/TEVC.2019.2938557.
%------------------------------- Copyright --------------------------------
% Copyright 2018-2019 Yiping Liu
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [-1,-1];
            obj.upper    = [1,1];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,~] = size(X);
            a     = 0.4;
            NP    = 2;
            psize = 0.10;
            center = [-0.50,0.50];
            Points = [center - psize; center + psize];
            PopObj = NaN(N,obj.M);
            for i = 1 : obj.M
                temp = abs(repmat(X(:,1),1,NP) - repmat(Points(i,:),N,1));
                temp(:,1) = temp(:,1) + 100 * (1 - cos(2*pi*(X(:,2)+0.5)));
                temp(:,2) = temp(:,2) + 100 * ((X(:,2)-0.5 + a*(X(:,1)-0.5)).^2);
                PopObj(:,i) = min(temp,[],2);
            end
            index = X(:,1) > 0;
            PopObj(index,:) = PopObj(index,:) + 0.01;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            div = 1/100;
            temp = 0:div:0.2;
            R = [temp',0.2-temp'];
        end
    end
end
