function P2 = Polygon_p_SampleInside(poly,n)
%Polygon_p_SampleInside Deterministic grid samples inside one polygon.
%   P2 = Polygon_p_SampleInside(poly,n) returns up to n points (m x 2) taken
%   from a uniform grid over the bounding box of poly ((M,2) vertices) and
%   kept if they are inside or on the boundary. If more than n points are
%   found they are subsampled deterministically by linspace indexing.

    px = poly(:,1);
    py = poly(:,2);
    xmin = min(px); xmax = max(px);
    ymin = min(py); ymax = max(py);
    ng   = max(10,ceil(3*sqrt(n)));
    [X,Y] = ndgrid(linspace(xmin,xmax,ng),linspace(ymin,ymax,ng));
    in   = inpolygon(X(:),Y(:),px,py);
    P2   = [X(in),Y(in)];
    if isempty(P2)
        P2 = [mean(px),mean(py)];
    end
    if size(P2,1) > n
        % unique() on the index list: round(linspace(...)) repeats indices
        % whenever the step is below 2, which would return duplicate rows and
        % collapse any nearest-neighbour based resolution estimate to zero.
        idx = unique(round(linspace(1,size(P2,1),n)));
        P2  = P2(idx,:);
    end
end
