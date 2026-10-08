function cand = CoSEANoveltyFilter(cand, Population, Problem, cfg)
%CoSEANoveltyFilter Keep only candidates far enough from the population.
%   cand = CoSEANoveltyFilter(cand, Population, Problem, cfg) drops injected
%   candidates whose minimum box-normalized decision-space distance to the
%   current population is below injEpsFrac*sqrt(d) (optional, default off).

    if ~cfg.noveltyFilter || isempty(cand)
        return;
    end
    lb = Problem.lower;
    ub = Problem.upper;
    sp = max(ub - lb, eps);
    Xn = Population.decs;                 % SOLUTION.decs takes no index
    Xn = (Xn - lb) ./ sp;
    Cn = (cand - lb) ./ sp;
    d  = Problem.D;
    Dm = pdist2(Cn, Xn);
    keep = min(Dm, [], 2).' > cfg.injEpsFrac * sqrt(d);
    cand = cand(keep, :);
end
