classdef Polygon_p < PROBLEM
% <2019> <multi> <real> <multimodal>
% Parameterized multi-modal multi-objective polygon distance minimization.
%   Parameters ('parameter',{M,D,row,col}):
%     M   - number of vertices of each polygon = number of objectives
%           (CalObj returns one column per vertex slot, minimized over the
%           row*col polygons elementwise, as in the original problem)
%     D   - number of decision variables (must be even; the 2-D layout is
%           tiled as repmat([x y],1,D/2), exactly as in the original)
%     row - number of polygons in a row    (default 2)
%     col - number of polygons in a column (default 2)
%   Number of equivalent Pareto sets (modes) = row*col.
%
%   Why this class exists: the platform port (Problems/.../Polygon/
%   PolygonProblem.m) hard-codes M=6, D=4, lower/upper=+-50 inside Setting,
%   and PROBLEM assigns external 'M'/'D' BEFORE Setting is called, so an
%   external sweep of M/D is silently overwritten. The original function
%   style implementation (Ref_algs_problems/PolygonTestProblem) uses
%   lower=-100, upper=100 and sweeps M in {3,4,6} and D in {4,8,10,...}.
%   This class restores the original bounds and makes M/D/row/col sweepable,
%   and adds GetPSModes (exact per-polygon Pareto sets) for mode-wise
%   recovery metrics.

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong. A multi-modal multi-objective
% evolutionary algorithm using two-archive and recombination strategies.
% IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
%--------------------------------------------------------------------------

    properties
        row;        % number of polygons in a row
        col;        % number of polygons in a column
        distance;   % distance between adjacent polygon centers
        Polygons;   % (M,2,row*col) vertices of each polygon
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [M,D,row,col] = obj.ParameterSet(6,4,2,2);
            obj.M        = M;
            obj.D        = D;
            obj.row      = row;
            obj.col      = col;
            obj.distance = 5;
            if mod(obj.D,2) ~= 0
                error('Polygon_p:oddD','D must be even (2-D layout is tiled).');
            end
            obj.lower    = -100*ones(1,obj.D);
            obj.upper    =  100*ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            obj.Polygons = Polygon_p_CreatePolygons(obj.row,obj.col,obj.distance,M);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            PopObj = Polygon_p_Obj(X,obj.Polygons,obj.D);
        end
        %% True Pareto sets per mode (one polygon = one mode)
        function modes = GetPSModes(obj,N)
            K      = obj.row * obj.col;
            nPer   = max(20,ceil(N/K));
            modes  = cell(1,K);
            for j = 1 : K
                P2 = Polygon_p_SampleInside(obj.Polygons(:,:,j),nPer);
                modes{j} = Polygon_p_Tile(P2,obj.D);
            end
        end
        %% Generate Pareto optimal solutions (objective space, PlatEMO convention)
        %  Sampled per mode (polygon) rather than over the whole decision box:
        %  the polygons cover only ~5% of the [-100,100]^2 layout extent, so a
        %  box grid would leave the reference PF far too sparse for IGD.
        function R = GetOptimum(obj,N)
            modes = obj.GetPSModes(N);
            PS    = vertcat(modes{:});
            R     = unique(Polygon_p_Obj(PS,obj.Polygons,obj.D),'rows');
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            modes = obj.GetPSModes(20000);
            PS    = vertcat(modes{:});
            R     = unique(Polygon_p_Obj(PS,obj.Polygons,obj.D),'rows');
        end
    end
end
