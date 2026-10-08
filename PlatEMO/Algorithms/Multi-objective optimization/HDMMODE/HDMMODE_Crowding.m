function CrowdDis = HDMMODE_Crowding(Pop,fk)
%HDMMODE_Crowding Harmonic average distance in the (normalized) decision space.
%   Ported from the original HDMMODE code (renamed with the HDMMODE_
%   prefix); logic unchanged.

    [N,~] = size(Pop);
    K = N - 1;
    Z    = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    pop = (Pop - repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance = pdist2(pop,pop,'minkowski',fk);
    value = sort(distance,2);
    CrowdDis = K./sum(1./value(:,2:N),2);
end
