function [PS,PF] = Polygon_p_PSnPF(Polygons,D,N)
%Polygon_p_PSnPF Reference Pareto set and front of the polygon problem.
%   [PS,PF] = Polygon_p_PSnPF(Polygons,D,N) samples a uniform 2-D grid over
%   the extent of all polygons, keeps the points inside any polygon (the true
%   Pareto set of the 2-D layout), tiles them into the D-dimensional decision
%   space with repmat([x y],1,D/2), and maps them to the objective space with
%   Polygon_p_Obj (the very same routine used by CalObj, so PS and PF are
%   consistent by construction).
%
%   Polygons - (M,2,K) vertices (output of Polygon_p_CreatePolygons)
%   D        - number of decision variables (even)
%   N        - target number of grid cells (>= number of returned points)

    [M,~,K] = size(Polygons); %#ok<ASGLU>
    xs = Polygons(:,1,:);
    ys = Polygons(:,2,:);
    xmin = min(xs(:))-1; xmax = max(xs(:))+1;
    ymin = min(ys(:))-1; ymax = max(ys(:))+1;

    ng = max(20,ceil(sqrt(N)));
    [X,Y] = ndgrid(linspace(xmin,xmax,ng),linspace(ymin,ymax,ng));
    X = X(:); Y = Y(:);
    inside = false(numel(X),1);
    for j = 1 : K
        inside = inside | inpolygon(X,Y,Polygons(:,1,j),Polygons(:,2,j));
    end
    P2 = [X(inside),Y(inside)];
    if size(P2,1) > N
        idx = unique(round(linspace(1,size(P2,1),N)));   % no duplicate rows
        P2  = P2(idx,:);
    end

    PS = Polygon_p_Tile(P2,D);
    PF = Polygon_p_Obj(PS,Polygons,D);
end
