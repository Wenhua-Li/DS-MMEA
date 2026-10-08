function s = CoSEASurvivalRate(candDecs, popDecs)
%CoSEASurvivalRate Fraction of injected candidates surviving the selection.
%   s = CoSEASurvivalRate(candDecs, popDecs) matches injected candidates to
%   the new population by near-exact decision vectors (injected candidates
%   pass through variation unchanged, so exact matching is valid) (M3).

    if isempty(candDecs) || isempty(popDecs)
        s = 0;
        return;
    end
    Dm  = pdist2(double(candDecs), double(popDecs));  % enforce double (single leaks from .mat)
    hit = min(Dm, [], 2) < 1e-10;
    s   = mean(hit);
end
