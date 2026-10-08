function PopObj = Polygon_p_Obj(X,Polygons,D)
%Polygon_p_Obj Objectives of the polygon problem for arbitrary X.
%   PopObj = Polygon_p_Obj(X,Polygons,D) computes, for each vertex slot k,
%     f_k(x) = min over polygons j of || x - tile(v_{j,k}) ||
%   where tile(v) = repmat(v,1,D/2). The distance is expanded analytically
%     ||x - tile(v)||^2 = Sxx + Syy - 2*(vx*Sx + vy*Sy) + C*(vx^2+vy^2)
%   with C = D/2 and Sx/Sy/Sxx/Syy the sums over the D/2 tiled copies of x,
%   so the result is numerically identical to the original tiled pdist2 call
%   for ANY x (not only for points on the tiled manifold), while costing
%   O(N*M*K) instead of building N x (M*C) distance matrices.
%
%   X        - n x D decision vectors
%   Polygons - (M,2,K) vertices (output of Polygon_p_CreatePolygons)
%   D        - number of decision variables (even)

    [M,~,K] = size(Polygons);
    C  = D/2;
    n  = size(X,1);
    Xr = reshape(X,n,2,C);
    Sx  = sum(Xr(:,1,:),3);
    Sy  = sum(Xr(:,2,:),3);
    Sxx = sum(Xr(:,1,:).^2,3);
    Syy = sum(Xr(:,2,:).^2,3);
    base = Sxx + Syy;                       % n x 1

    PopObj = Inf(n,M);
    for j = 1 : K
        vx = reshape(Polygons(:,1,j),1,M);  % 1 x M
        vy = reshape(Polygons(:,2,j),1,M);
        d2 = base - 2*(Sx.*vx + Sy.*vy) + C*(vx.^2 + vy.^2);
        d2(d2 < 0) = 0;                     % guard round-off
        PopObj = min(PopObj,sqrt(d2));
    end
end
