function Crowd = FPITSEA_Crowding(Pop)
% Harmonic average distance of each solution in the decision space.
% Renamed from Crowding to avoid potential conflicts across algorithms.

    [N,~] = size(Pop);
    K     = N - 1;
    Z     = min(Pop,[],1);
    Zmax  = max(Pop,[],1);
    Norpop = (Pop - repmat(Z,N,1)) ./ repmat(Zmax - Z,N,1);
    distance = pdist2(Norpop,Norpop);
    [value,~] = sort(distance);
    Crowd = K ./ sum(1./value(2:N,:));
end