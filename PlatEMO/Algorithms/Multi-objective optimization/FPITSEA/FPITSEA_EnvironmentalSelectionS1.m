function [Population,value] = FPITSEA_EnvironmentalSelectionS1(Union,N)
% Environmental selection in the first stage of FPITSEA.

    %% Calculate fuzzy preference indicator
    IFP = FPITSEA_CalFP(Union);

    %% Select the next generation based on the fuzzy preference indicator
    [~,idx]      = sort(IFP,'descend');
    Population   = Union(idx(1:N));

    %% Calculate crowding distance for mating selection
    value = CalCD(Population.decs);
end

function CD = CalCD(Pop)
% Crowding distance used for mating selection.
    N    = size(Pop,1);
    Dist = pdist2(Pop,Pop);
    Dist(logical(eye(size(Pop,1)))) = inf;
    Dist = sort(Dist,2);
    CD   = 1 ./ (sum(Dist(:,1:floor(sqrt(N))),2) + 2);
end