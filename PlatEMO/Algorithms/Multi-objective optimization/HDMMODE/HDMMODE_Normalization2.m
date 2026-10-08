function outcome = HDMMODE_Normalization2(Population,Average,G)
%HDMMODE_Normalization2 Record the mean of the min-max normalized objectives.
%   Ported from the original HDMMODE code (renamed); logic unchanged.

    fmax = max(Population.objs,[],1);
    fmin = min(Population.objs,[],1);
    Average(G,:) = mean((Population.objs - fmin)./(fmax - fmin));
    outcome = Average;
end
