function MatingPool = CoSEAMateSelect(Population, FrontNo, Div, Problem, cfg)
%CoSEAMateSelect Mating selection of CoSEA.
%   Default: binary tournament on (FrontNo, -Div) as in NSGA-II, with the
%   decision-space diversity measure Div replacing crowding distance.
%   Optional (cfg.mateRestrict, M8): each candidate first competes with its
%   decision-space nearest neighbor (DN-NSGA-II style implicit intra-niche
%   mating restriction), the winner enters the tournament bracket.

    N = Problem.N;
    if ~cfg.mateRestrict
        MatingPool = TournamentSelection(2, N, FrontNo, -Div);
        return;
    end
    n  = length(Population);
    Xn = Population.decs;                 % SOLUTION.decs takes no index
    lb = Problem.lower;
    ub = Problem.upper;
    Xn = (Xn - lb) ./ max(ub - lb, eps);
    D  = pdist2(Xn, Xn);
    D(1:n+1:end) = inf;
    [~, nnIdx] = min(D, [], 2);
    Pool = zeros(1, N);
    for i = 1 : N
        a = randi(n);
        b = nnIdx(a);
        if FrontNo(a) < FrontNo(b) || (FrontNo(a) == FrontNo(b) && Div(a) >= Div(b))
            Pool(i) = a;
        else
            Pool(i) = b;
        end
    end
    MatingPool = Pool;
end
