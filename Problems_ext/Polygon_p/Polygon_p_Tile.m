function XD = Polygon_p_Tile(P2,D)
%Polygon_p_Tile Embed 2-D points into the D-dimensional decision space.
%   XD = Polygon_p_Tile(P2,D) maps each row [x y] of P2 (n x 2) to
%   repmat([x y],1,D/2), the tiling used by the original polygon problem
%   (transform.m). D must be even.

    if mod(D,2) ~= 0
        error('Polygon_p_Tile:oddD','D must be even.');
    end
    XD = repmat(P2,1,D/2);
end
