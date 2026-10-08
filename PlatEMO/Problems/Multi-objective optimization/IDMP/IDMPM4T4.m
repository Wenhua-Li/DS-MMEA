classdef IDMPM4T4 < PROBLEM
% <2019> <multi> <real> <multimodal>
% Imbalanced distance minimization problem with 4 objectives and 4 segments

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, G. G. Yen, Y. Nojima and N. Masuyama, Handling
% Imbalance Between Convergence and Diversity in the Decision Space in
% Evolutionary Multimodal Multi-Objective Optimization, IEEE Transactions
% on Evolutionary Computation, 2019, Early Access,
% DOI: 10.1109/TEVC.2019.2938557.
%------------------------------- Copyright --------------------------------
% Copyright 2018-2019 Yiping Liu (port to PlatEMO 2026 class architecture)
%--------------------------------------------------------------------------

    properties(Access = private)
        Points;
    end
    methods
        function Setting(obj)
            obj.M        = 4;
            obj.D        = 4;
            obj.lower    = [-1,-1,-1,-1];
            obj.upper    = [1,1,1,1];
            obj.encoding = ones(1,obj.D);
            pgon    = nsidedpoly(obj.M);
            psize   = [.1,.1,.1,.1];
            center  = [-.50,-.50;.50,-.50;.50,.50;-.50,.50];
            obj.Points = NaN(obj.M,2,4);
            for i = 1 : 4
                obj.Points(:,:,i) = pgon.Vertices.*psize(i) + center(i,:);
            end
        end
        function PopObj = CalObj(obj,X)
            [N,~] = size(X);
            PopObj = NaN(N,obj.M);
            for i = 1 : obj.M
                temp = pdist2(X(:,1:2),reshape(obj.Points(i,:,:),[2,4])');
                temp(:,1) = temp(:,1) + 100*((X(:,3)+.6).^2+1-cos(1*2*pi*(X(:,3)+.6))) + 100*((X(:,4)+.6).^2+1-cos(0*2*pi*(X(:,4)+.6)));
                temp(:,2) = temp(:,2) + 100*((X(:,3)+.2).^2+1-cos(2*2*pi*(X(:,3)+.2))) + 100*((X(:,4)+.2).^2+1-cos(0*2*pi*(X(:,4)+.2)));
                temp(:,3) = temp(:,3) + 100*((X(:,3)-.2).^2+1-cos(3*2*pi*(X(:,3)-.2))) + 100*((X(:,4)-.2).^2+1-cos(0*2*pi*(X(:,4)-.2)));
                temp(:,4) = temp(:,4) + 100*((X(:,3)-.6).^2+1-cos(4*2*pi*(X(:,3)-.6))) + 100*((X(:,4)-.6).^2+1-cos(0*2*pi*(X(:,4)-.6)));
                PopObj(:,i) = min(temp,[],2);
            end
        end
        function R = GetPF(obj)
            [X,Y] = ndgrid(linspace(0,1,100));
            ND    = inpolygon(X(:),Y(:),obj.Points(:,1,3),obj.Points(:,2,3));
            R     = pdist2([X(ND),Y(ND)],obj.Points(:,:,3));
        end
    end
end