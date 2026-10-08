function Polygons = CreatePolygons(row,col,distance,M)
% Create 2D regular polygons arranged in a grid.
%
% row      - Number of polygons in a row
% col      - Number of polygons in a column
% distance - Distance between adjacent polygon centers
% M        - Number of vertices (equals the number of objectives)
%
% Polygons - (M,2,N) array, where Polygons(:,:,i) gives the vertices of
%            the i-th polygon.

    Centers = zeros(row*col,2);
    cnt = 1;
    for i = 1 : row
        for j = 1 : col
            Centers(cnt,:) = [distance*(j-1), distance*(i-1)];
            cnt = cnt + 1;
        end
    end

    scenario = 1;
    if scenario == 1
        Angle = 2*pi.*(1:M)/M;
    else
        if mod(M,2) == 0
            Angle = (2.*(1:M)-3).*pi./M;
        else
            Angle = (2.*(1:M)-2).*pi./M;
        end
    end

    N = row*col;
    Polygons = zeros(length(Angle),2,N);
    for i = 1 : N
        vertices = ([sin(Angle)', cos(Angle)']) + repmat(Centers(i,:),length(Angle),1);
        Polygons(:,:,i) = vertices;
    end
end