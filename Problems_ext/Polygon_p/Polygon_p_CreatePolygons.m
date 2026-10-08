function Polygons = Polygon_p_CreatePolygons(row,col,distance,M)
%Polygon_p_CreatePolygons Regular M-gons arranged in a row x col grid.
%   Polygons = Polygon_p_CreatePolygons(row,col,distance,M) returns an
%   (M,2,row*col) array; Polygons(:,:,j) holds the M vertices of polygon j.
%   Vertex angles and center layout are identical to the original
%   CreatePolygons.m (scenario = 1), only the function name is unique to
%   avoid same-name resolution problems under addpath(genpath(...)).

    Centers = zeros(row*col,2);
    cnt = 1;
    for i = 1 : row
        for j = 1 : col
            Centers(cnt,:) = [distance*(j-1), distance*(i-1)];
            cnt = cnt + 1;
        end
    end

    Angle = 2*pi.*(1:M)/M;

    K = row*col;
    Polygons = zeros(M,2,K);
    for i = 1 : K
        Polygons(:,:,i) = [sin(Angle)', cos(Angle)'] + repmat(Centers(i,:),M,1);
    end
end
