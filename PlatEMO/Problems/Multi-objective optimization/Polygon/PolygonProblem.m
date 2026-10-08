classdef PolygonProblem < PROBLEM
% <2018> <multi> <real> <multimodal>
% The multi-modal multi-objective polygon distance minimization problem
% lower    --- -50 --- Lower bound of decision variables
% upper    ---  50 --- Upper bound of decision variables
% row      ---   2  --- Number of polygons in a row
% col      ---   2  --- Number of polygons in a column
% distance ---   5  --- Distance between adjacent polygon centers

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong. A multi-modal multi-objective
% evolutionary algorithm using two-archive and recombination strategies.
% IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
%--------------------------------------------------------------------------

    properties
        row;       % Number of polygons in a row
        col;       % Number of polygons in a column
        distance;  % Distance between polygon centers
        Polygons;  % Vertices of each polygon (M,2,Np)
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.row,obj.col,obj.distance] = obj.ParameterSet(2,2,5);
            obj.M        = 6;
            obj.D        = 4;
            obj.lower    = -50*ones(1,obj.D);
            obj.upper    =  50*ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.Polygons = CreatePolygons(obj.row,obj.col,obj.distance,obj.M);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            Np = obj.row * obj.col;
            PopObj = Inf;
            for i = 1 : Np
                % Map 2D polygon to D-dimensional decision space using two basis vectors
                vertices = repmat(obj.Polygons(:,:,i),1,obj.D/2);
                PopObj   = min(PopObj,pdist2(X,vertices));
            end
        end
        %% Sample points on the Pareto optimal set and the Pareto front
        function R = GetOptimum(obj,N)
            % For multimodal MOO, the optimum returned here is used as PF (objective values);
            % the Pareto-optimal set is reconstructed via PSnPF on demand.
            [~,PF] = PSnPF(obj.Polygons,obj.D,N);
            R = PF;
        end
        function R = GetPF(obj)
            [~,PF] = PSnPF(obj.Polygons,obj.D,100);
            R = PF;
        end
    end
end