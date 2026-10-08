function [Population, MRank, Div, info] = MRCoSEAEnvSelect(Population, N, Problem, cfg, progress)
%MRCoSEAEnvSelect Mode-relative environmental selection (the v3.0 operator).
%   [Population, MRank, Div, info] = MRCoSEAEnvSelect(...) implements the
%   single v3.0 selection principle: an individual's survival is decided by
%   its Pareto rank WITHIN ITS OWN MODE, subject to one annealed
%   epsilon-relaxation gate against the global front. This one operator
%   replaces the v2.0 stack (per-mode M8 epsilon, local front-1
%   reservation, persistent archive; scheduling and quotas lived in the
%   main loop and are deleted in v3.0).
%
%     groups   - the merged population is clustered by MRCoSEAModeCluster;
%                label >= 1 members are ranked by within-mode non-dominated
%                sorting (a local-PF mode is front-1 within its mode and is
%                never forced to compete against the global front), label 0
%                (unclustered) members keep their GLOBAL front number (a
%                reflected candidate landing on an unvisited Pareto set is
%                globally front-1 and protected before any mode exists).
%     MRank    - within-group front number of every individual.
%     gate     - group front-1 members are protected unless dominated by
%                the min-max-normalized global front inflated by epsT
%                (z_front + epsT <= z_i componentwise). epsT anneals from
%                epsAnneal0 to epsFloor (convergence pressure + junk-blob
%                pruning; a local-PF mode within epsFloor of the front
%                survives indefinitely).
%     overflow - if the protected set exceeds N slots, the most isolated
%                member of every mode is kept first (a small mode can never
%                be truncated away entirely), the rest is truncated by
%                SPEA2 in the decision space (decision-space coverage of
%                the protected set directly targets IGDX/CR uniformity).
%     layer 2  - remaining slots are filled by the v0.8 global-front path
%                (complete fronts first, SPEA2 truncation at the boundary).
%   info reports protN (protected count), K (mode count) and epsT.

    X    = Population.decs;               % SOLUTION.decs takes no index (E2)
    ObjM = Population.objs;               % same for .objs/.cons
    ConM = Population.cons;
    nAll = size(X, 1);
    lb   = Problem.lower;
    ub   = Problem.upper;
    Xn   = (X - lb) ./ max(ub - lb, eps); % box normalization (B3: no drift)

    %% Group ranking: within-mode fronts, noise by global fronts
    lab = MRCoSEAModeCluster(X, Problem, cfg);
    uL  = unique(lab(lab > 0));
    [FrontNo, ~] = NDSort(ObjM, ConM, N);
    FrontNo = FrontNo(:)';
    Fidx    = find(FrontNo == 1);
    epsT    = MRCoSEAEpsAnneal(progress, cfg);

    %% Objective normalization for the gate / mode-net. v3.0 stretches the
    %  scale by the merged population min/max; objective-space outliers
    %  (NMMF2 f2 up to ~1e5, MMF1_e) stretch it 30-1000x and neutralize the
    %  gate (Stage-0 diagnosis). cfg.robustNorm bases the scale on the
    %  global front-1 range instead. The default keeps v3.0 semantics.
    if cfg.robustNorm && ~isempty(Fidx)
        zmin = min(ObjM(Fidx,:), [], 1);
        zmax = max(ObjM(Fidx,:), [], 1);
    else
        zmin = min(ObjM, [], 1);
        zmax = max(ObjM, [], 1);
    end
    Zn = (ObjM - zmin) ./ max(zmax - zmin, eps);

    MRank = FrontNo;                    % default: global (label 0)
    useModeNet = strcmpi(cfg.gateType, 'modeNet');
    for k = 1 : numel(uL)
        idx = find(lab == uL(k));
        if useModeNet
            % v3.1 arm: within-mode epsilon-grid fronts (grid surrogate of
            % additive eps-dominance, as in eps-MOEA). Near-tie compression
            % inside each mode liberates overflow slots for small modes;
            % junk modes keep their own front (accepted by design). The
            % grid scale is 0.1x the global anneal clamped to [0.02, 0.06]
            % (documented constant): within-mode gaps are far smaller than
            % the cross-mode margins the global curve was designed for.
            epsM = min(max(0.1 * epsT, 0.02), 0.06);
            Zg   = floor(Zn(idx,:) ./ max(epsM, realmin));
            fno  = NDSort(Zg, numel(idx));
        else
            fno = NDSort(ObjM(idx,:), ConM(idx,:), numel(idx));
        end
        MRank(idx) = fno(:)';             % within-mode fronts
    end
    switch cfg.envSel
        case 'crowding'
            Div = CoSEACrowdingDec(Xn, MRank);
        otherwise
            Div = CoSEAKnnIsolation(Xn, MRank);
    end

    %% Protected set: group front-1 members; junk pruning per cfg.gateType
    Prot   = (MRank == 1);
    gated  = 0;
    floorSkip = 0;
    floorAdd  = 0;
    switch lower(cfg.gateType)
        case {'modenet', 'none'}
            % no global pruning: 'modenet' additionally ranks members by
            % the within-mode eps grid (useModeNet above); 'none' keeps the
            % exact within-mode fronts (the pure mode-relative arm, v3.1)
        case 'dualref'
            % v3.1 arm: mode-level global judgment. A mode is pruned only
            % WHOLESALE: if even its ideal point (componentwise best over
            % all members) is eps-dominated by the global front, none of
            % its front-1 members survives. A mode anchored by any global
            % front-1 member is exempt, so a developing colony holding one
            % converged member can no longer be picked off member by member
            % (the 17x/11x disasters on IDMPM4T1/T3 in the blame map).
            for k = 1 : numel(uL)
                mk  = find(lab == uL(k));
                mk1 = mk(MRank(mk) == 1);
                if isempty(mk1) || any(FrontNo(mk1) == 1)
                    continue;
                end
                ideal = min(Zn(mk,:), [], 1);
                if any(all(Zn(Fidx,:) + epsT <= ideal, 2))
                    Prot(mk1) = false;
                    gated     = gated + numel(mk1);
                end
            end
        otherwise  % 'global': the v3.0 per-member gate (bit-identical)
            if any(Prot)
                cand = find(Prot & FrontNo > 1);  % global front passes free
                for t = 1 : numel(cand)
                    c = cand(t);
                    if any(all(Zn(Fidx,:) + epsT <= Zn(c,:), 2))
                        Prot(c) = false;
                        gated   = gated + 1;
                    end
                end
            end
    end

    %% Ungated per-mode floor (cfg.modeFloor): every mode keeps its most
    % isolated front-1 representatives regardless of the gate. A developing
    % mode whose whole front-1 is still outside the band must not lose every
    % member; the v2.4/v2.5 experiment on the v2.0 reservation showed that
    % gating developing modes backfires (they are indistinguishable from
    % local-PF modes in the delta criterion), so the floor is unconditional.
    % v3.0 had deleted the v2.0 unconditional reservation together with the
    % scheduling stack; without the floor the weight-imbalanced weak modes
    % of IDMPT1 die (G5 regression, pilot P3.1 diagnosis).
    if cfg.modeFloor > 0
        for k = 1 : numel(uL)
            mk = find(lab == uL(k));
            % Size guard: the floor exists for DEVELOPING small modes; a
            % mode that holds a large share of the population (the
            % degenerate single-blob regime of high-dimensional Polygon)
            % does not need it and protecting it stalls convergence.
            if numel(mk) > cfg.modeFloorMaxFrac * N
                floorSkip = floorSkip + 1;
                continue;
            end
            mk = mk(MRank(mk) == 1 & ~Prot(mk));
            if isempty(mk)
                continue;
            end
            [~, ord] = sort(Div(mk), 'descend');
            take = mk(ord(1:min(cfg.modeFloor, numel(mk))));
            Prot(take) = true;
            floorAdd = floorAdd + numel(take);
        end
    end

    %% Assemble: protected first (per-mode guarantee + SPEA2), then layer 2
    Next = false(1, nAll);
    nP   = sum(Prot);
    if nP > 0
        if nP <= N
            Next(Prot) = true;
        else
            % Keep the most isolated member of every mode first so that a
            % small mode can never be truncated away entirely.
            for k = 1 : numel(uL)
                mk = find(Prot & lab == uL(k));
                if isempty(mk)
                    continue;
                end
                [~, b] = max(Div(mk));
                Next(mk(b)) = true;
            end
            rest = find(Prot & ~Next);
            slotsP = N - sum(Next);
            if ~isempty(rest) && slotsP > 0
                keep = CoSEASpeaTrunc(Xn(rest,:), slotsP);
                Next(rest(keep)) = true;
            end
            % slotsP == 0 with rest non-empty: already full (more modes
            % than slots), nothing more to keep
        end
    end

    %% Layer 2: v0.8 global-front fill of the remaining slots
    slots = N - sum(Next);
    if slots > 0
        Rem = find(~Next);
        if ~isempty(Rem)
            [~, ord] = sort(FrontNo(Rem));
            Rem  = Rem(ord);
            fRem = FrontNo(Rem);
            uF   = unique(fRem);
            for u = 1 : numel(uF)
                idxU = Rem(fRem == uF(u));
                if numel(idxU) <= slots
                    Next(idxU) = true;
                    slots = slots - numel(idxU);
                else
                    if strcmp(cfg.envSel, 'crowding')
                        [~, Rank] = sort(Div(idxU), 'descend');
                        Next(idxU(Rank(1:slots))) = true;
                    else
                        keep = CoSEASpeaTrunc(Xn(idxU,:), slots);
                        Next(idxU(keep)) = true;
                    end
                    slots = 0;
                end
                if slots <= 0
                    break;
                end
            end
        end
    end

    %% Population for the next generation
    Population = Population(Next);
    MRank      = MRank(Next);
    Div        = Div(Next);
    maxModeSize = 0;
    for k = 1 : numel(uL)
        maxModeSize = max(maxModeSize, sum(lab == uL(k)));
    end
    info       = struct('protN', nP, 'K', numel(uL), 'epsT', epsT, ...
        'gated', gated, 'floorSkip', floorSkip, 'floorAdd', floorAdd, ...
        'maxModeSize', maxModeSize, ...
        'lab', lab(Next));  % survivor cluster labels (report layer reads
                            % these; adding the field is behavior-neutral)
end

function epsT = MRCoSEAEpsAnneal(progress, cfg)
%MRCoSEAEpsAnneal Gate epsilon schedule on the normalized objectives.
%   'linear' : epsFloor + (epsAnneal0 - epsFloor)*(1 - progress)
%   'commea' : CoMMEA-style -log2(1.5*t) clamped to [epsFloor, epsAnneal0]
%              (reaches the floor at t = 2^(-epsFloor)/1.5, e.g. 0.40 for
%              epsFloor = 0.25)
    switch lower(cfg.annealShape)
        case 'commea'
            epsT = -log2(1.5 * max(progress, 1e-12));
            epsT = min(max(epsT, cfg.epsFloor), cfg.epsAnneal0);
        otherwise  % 'linear'
            epsT = cfg.epsFloor + ...
                (cfg.epsAnneal0 - cfg.epsFloor) * max(0, 1 - progress);
    end
end
