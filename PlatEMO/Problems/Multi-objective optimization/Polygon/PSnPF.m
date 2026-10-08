function [PS,PF] = PSnPF(Polygons,D,N)
% Sample reference points for the polygon test problem.
%
% Polygons - (M,2,K) polygons (output of CreatePolygons)
% D        - Dimension of the decision space
% N        - Number of reference points to sample

    [M,~,K] = size(Polygons);

    % Obtain PS by sampling a 2D grid and keeping points inside any polygon
    [X,Y] = ndgrid(linspace(-2,10,ceil(sqrt(N))));
    nondominated = zeros(numel(X),1);
    for i = 1 : K
        polygon = Polygons(:,:,i);
        ND = inpolygon(X(:),Y(:),polygon(:,1),polygon(:,2));
        nondominated = nondominated | ND;
    end
    PS = [X(nondominated), Y(nondominated)];

    % Obtain PF by mapping each PS point to objective space
    num_PS = size(PS,1);
    objectives = Inf(num_PS,M);
    for i = 1 : num_PS
        point = PS(i,:);
        for j = 1 : K
            dist = pdist2(transform(point,D), transform(Polygons(:,:,j),D));
            objectives(i,:) = min(objectives(i,:), dist);
        end
    end
    PF = objectives;
end

function B = transform(A,D)
% Map a 2D point A to D-dimensional space by tiling with two basis vectors.
    B = repmat(A,1,D/2);
end