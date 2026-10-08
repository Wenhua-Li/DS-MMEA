function modes = CoSEAClusterModes(PS,lower,upper,thr)
%CoSEAClusterModes Split a true Pareto set sample into equivalent PSs (modes).
%   modes = CoSEAClusterModes(PS,lower,upper,thr) clusters the rows of PS by
%   connectivity in the box-normalized decision space: two points are linked
%   when their normalized Euclidean distance is below thr (default
%   0.05*sqrt(D), the convention of CoSEATruePSModes). Each connected
%   component is one mode. Modes are returned in increasing order of the
%   centroid of their first decision variable (deterministic).
%   For large samples (> 2500 points) the connectivity graph is built from a
%   k-nearest-neighbour search (k = 12) instead of the full distance matrix,
%   which is equivalent for the benchmark PS geometries used here (each mode
%   is a densely sampled connected manifold) and avoids O(n^2) memory.

    if nargin < 4 || isempty(thr)
        thr = 0.05*sqrt(size(PS,2));
    end
    n = size(PS,1);
    if n < 2
        modes = {PS};
        return;
    end

    sp = max(upper - lower,eps);
    Pn = (PS - lower) ./ sp;

    if n <= 2500
        Dm  = pdist2(Pn,Pn);
        adj = sparse(Dm < thr);
        adj = adj | adj';
    else
        k    = min(12,n-1);
        [idx,dist] = knnsearch(Pn,Pn,'K',k+1);
        rows = repmat((1:n)',1,k+1);
        keep = dist < thr;
        adj  = sparse(rows(keep),idx(keep),true,n,n);
        adj  = adj | adj';
    end
    adj(1:n+1:end) = false;                 % no self loops

    comp  = conncomp(graph(adj));
    nComp = max(comp);
    modes = cell(1,nComp);
    for c = 1 : nComp
        modes{c} = PS(comp == c,:);
    end

    cen = zeros(1,nComp);
    for c = 1 : nComp
        cen(c) = mean(modes{c}(:,1));
    end
    [~,ord] = sort(cen);
    modes = modes(ord);
end
