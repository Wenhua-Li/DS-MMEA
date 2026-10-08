function [Population, FrontNo, Div] = CoSEAEnvSelect(Population, N, Problem, cfg, progress)
%CoSEAEnvSelect Environmental selection of CoSEA (self-contained).
%   Convergence: fast non-dominated sorting in the objective space (NDSort).
%   Diversity:   decision-space maintenance, three paths:
%     'crowding' - v0.6 axis-aligned crowding distance (ablation only)
%     'spea'     - v0.7 kNN isolation (tournament) + SPEA2 truncation
%     'speaEps'  - v0.8 = v0.7 + local-PS retention (epsilon-dominance, M8)
%   Returns the selected population, front numbers and the diversity measure
%   (isolation or crowding, larger = more isolated = preferred in tournament).

    [FrontNo, MaxFNo] = NDSort(Population.objs, Population.cons, N);
    X  = Population.decs;                 % SOLUTION.decs takes no index (E2)
    lb = Problem.lower;
    ub = Problem.upper;
    Xn = (X - lb) ./ max(ub - lb, eps);   % box normalization (B3: no drift)

    switch cfg.envSel
        case 'crowding'
            Div = CoSEACrowdingDec(Xn, FrontNo);
        otherwise % 'spea' | 'speaEps'
            Div = CoSEAKnnIsolation(Xn, FrontNo);
    end

    Next = FrontNo < MaxFNo;

    %% M8: local-PS retention (epsilon-dominance), decays over the budget
    retain = false(size(FrontNo));
    if strcmp(cfg.envSel, 'speaEps')
        epsT   = cfg.epsRetainMin + (cfg.epsRetain0 - cfg.epsRetainMin) * (1 - progress);
        retain = CoSEAEpsRetain(Population.objs, FrontNo, epsT);
        retain(Next) = false;             % already selected by rank
        cap = ceil(cfg.epsRetainCapFrac * N);
        if sum(retain) > cap              % keep the most isolated retained
            idxR = find(retain);
            [~, ord] = sort(Div(idxR), 'descend');
            keepR  = idxR(ord(1:cap));
            retain = false(size(FrontNo));
            retain(keepR) = true;
        end
    end

    %% Fill the remaining slots from the critical front
    slots = N - sum(Next) - sum(retain);
    if slots < 0
        % Protected fronts plus retained exceed N: drop least isolated retained
        idxR = find(retain);
        [~, ord] = sort(Div(idxR), 'ascend');
        retain(idxR(ord(1:-slots))) = false;
        slots = 0;
    end
    Last = find(FrontNo == MaxFNo & ~retain);
    if ~isempty(Last) && slots > 0
        if strcmp(cfg.envSel, 'crowding')
            [~, Rank] = sort(Div(Last), 'descend');
            chosen = Last(Rank(1:min(slots, numel(Last))));
        else
            nTake = min(slots, numel(Last));
            if nTake < numel(Last)
                keepL = CoSEASpeaTrunc(Xn(Last,:), nTake);
                chosen = Last(keepL);
            else
                chosen = Last;
            end
        end
        Next(chosen) = true;
    end
    Next(retain) = true;

    %% Population for the next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    Div        = Div(Next);
end
