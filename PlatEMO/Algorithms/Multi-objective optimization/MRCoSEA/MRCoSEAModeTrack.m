function [Modes, nextId] = MRCoSEAModeTrack(Modes, nextId, lab, Xn, cfg)
%MRCoSEAModeTrack Cross-generation mode tracking (id, centroid, size, age).
%   [Modes, nextId] = MRCoSEAModeTrack(...) matches the current mode
%   clusters (labels lab from MRCoSEAModeCluster, ordered by decreasing
%   size) to the previous mode list by nearest centroid (greedy,
%   deterministic, radius matchThrFrac*epsFrac*sqrt(d)). Matched modes
%   inherit their persistent id and age (+1); unmatched clusters are
%   newborn (age 0, fresh id). v3.0 keeps ONLY this lightweight tracking
%   (for the maturity gate of the manifold channel and the per-mode
%   survival feedback); the v2.0 indicator machinery (three-state machine,
%   LocalC, eps-progress) is deleted.

    d   = size(Xn, 2);
    thr = cfg.matchThrFrac * cfg.epsFrac * sqrt(d);

    %% Current clusters (labels 1..K, size-ordered by MRCoSEAModeCluster)
    K    = max([lab, 0]);
    cent = zeros(K, d);
    sz   = zeros(1, K);
    for k = 1 : K
        idx = lab == k;
        cent(k, :) = mean(Xn(idx, :), 1);
        sz(k) = sum(idx);
    end

    %% Greedy nearest-centroid matching old -> new (deterministic)
    match = zeros(1, K);                          % new k <- old match(k)
    if ~isempty(Modes) && K > 0
        C = pdist2(vertcat(Modes.cent), cent);
        usedOld = false(1, numel(Modes));
        while true
            [cmin, idx] = min(C(:));
            if cmin >= thr
                break;
            end
            [jo, jn] = ind2sub(size(C), idx);
            if usedOld(jo) || match(jn) ~= 0
                C(jo, jn) = inf;
                continue;
            end
            match(jn) = jo;
            usedOld(jo) = true;
            C(jo, :) = inf;
            C(:, jn) = inf;
        end
    end

    %% Build the new mode list (field order fixed for struct arrays)
    newModes = struct('id',{},'cent',{},'size',{},'age',{}); %#ok<STRUCT>
    for k = 1 : K
        j = match(k);
        if j == 0                                   % newborn mode
            nextId = nextId + 1;
            mk = struct('id',nextId,'cent',cent(k,:),'size',sz(k),'age',0);
        else
            mk = struct('id',Modes(j).id,'cent',cent(k,:),'size',sz(k), ...
                'age',Modes(j).age + 1);
        end
        newModes(end+1) = mk; %#ok<AGROW>
    end
    Modes = newModes;
end
