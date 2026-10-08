classdef MRCoSEA < ALGORITHM
% <2026> <multi> <real> <multimodal>
% MRCoSEA - Mode-Relative CoSEA (CoSEA v3.0) for multimodal multiobjective
%   optimization. One selection principle replaces the five v2.0 patch
%   mechanisms (per-mode epsilon, local front-1 reservation, persistent
%   archive, indicator scheduling, offspring quotas):
%     survival - an individual's survival is decided by its Pareto rank
%                WITHIN ITS OWN MODE (MRCoSEAEnvSelect): mode front-1
%                members survive structurally (a local-PF mode is never
%                forced to compete against the global front), subject to a
%                single annealed epsilon-relaxation gate against the global
%                front that provides convergence pressure and prunes junk
%                blobs. Unclustered points (label 0) are ranked by their
%                global front number, so a reflected candidate landing on
%                an unvisited Pareto set is globally front-1 and protected
%                before any mode exists there (discovery channel).
%     discover - manifold-aligned symmetric reflection injection with M3
%                survival-rate feedback and auto-degradation (unchanged
%                from CoSEA v0.8).
%     exploit  - per-mode local-PCA manifold sampling (RM-MEDA style,
%                closed form) with per-mode survival-rate feedback, on the
%                unchanged v0.8 backbone (non-dominated sorting is only
%                used inside the mode-relative operator and the fill path).
%   Mating is a plain global binary tournament on (mode rank, -isolation):
%   with survival safe at the selection level, restrictive mating is no
%   longer needed (v2.0 pool mating caused the G1 mixing loss).
%   All offspring and injected candidates are really evaluated and counted
%   in the budget (Problem.Evaluation), no free evaluations.
%   With mrEnabled=false and manifoldEnabled=false the algorithm is
%   bit-identical to CoSEA v0.8 (the off-path calls the v0.8 functions
%   directly; regression-gated by exp/mr/test_mr_equiv.m).
%
%   Configuration: MRCoSEA('parameter',{struct('mrEnabled',false,...)})
%   overrides any subset of MRCoSEAConfig defaults (no workspace globals,
%   parfor-safe).

%------------------------------- Copyright --------------------------------
% Thesis project 2026TEVC S2 (CoSEA v3.0 mode-relative revision, 2026-09).
%--------------------------------------------------------------------------

    properties
        ckpt = {};  % checkpoints: {FE, PopDec, PopObj} at every 20% budget
        log  = struct('inj',[],'surv',[],'m',[], ...      % v0.8 injection log
            'gen',[],'nModes',[],'maniN',[],'protN',[], ... % per-generation MR
            'maxModeSize',[],'floorSkip',[],'floorAdd',[], ... % floor diagnostics
            'modeTab',zeros(0,5)); % [gen id size age qModel]
    end
    methods
        function main(Algorithm, Problem)
            cfg = MRCoSEAConfig(Algorithm.ParameterSet(struct()));
            N     = Problem.N;
            maxFE = Problem.maxFE;

            %% MR switches: any new feature on reroutes the main loop
            mrOn  = cfg.mrEnabled;
            newOn = cfg.mrEnabled || cfg.manifoldEnabled;

            %% Scheduling constants (in generations, v0.8 semantics)
            maxgen   = ceil((maxFE - N) / N);
            kInj     = max(2, round(cfg.kFrac * maxgen));
            startGen = ceil(cfg.injectStartFrac * maxgen);
            if isfield(cfg,'injectStopFrac') && cfg.injectStopFrac < 1
                stopGen = floor(cfg.injectStopFrac * maxgen);
            else
                stopGen = Inf;
            end

            %% M3 state (v0.8)
            mMax       = round(cfg.mFrac * N);
            mMin       = max(1, round(cfg.mMinFrac * N));
            mAdapt     = mMax;
            injEnabled = true;
            zeroCnt    = 0;
            lastCand   = [];

            %% Mode state (MR): persistent mode list for age tracking and
            %% the per-mode model-share feedback (no indicators, no classes)
            Modes  = struct('id',{},'cent',{},'size',{},'age',{}); %#ok<STRUCT>
            fracA  = struct('id',{},'f',{});          % adaptive model shares
            nextId = 0;

            %% Generate random population
            Population = Problem.Initialization();
            if newOn
                if mrOn
                    [Population, MRank, Div] = MRCoSEAEnvSelect( ...
                        Population, N, Problem, cfg, 0);
                else
                    [Population, FrontNo, Div] = CoSEAEnvSelect( ...
                        Population, N, Problem, cfg, 0);
                end
            else
                [Population, FrontNo, Div] = CoSEAEnvSelect( ...
                    Population, N, Problem, cfg, 0);
            end
            gen      = 0;
            nextCkpt = 0.2 * maxFE;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                progress = Problem.FE / maxFE;

                if ~newOn
                    %% v0.8 path (bit-identical fallback)
                    MatingPool = CoSEAMateSelect(Population, FrontNo, Div, Problem, cfg);
                    Offspring  = OperatorGA(Problem, Population(MatingPool));
                else
                    %% MR path: track modes on the parents
                    X   = Population.decs;        % E2: no indexed .decs
                    lab = MRCoSEAModeCluster(X, Problem, cfg);
                    Xn  = (X - Problem.lower) ./ ...
                          max(Problem.upper - Problem.lower, eps);
                    [Modes, nextId] = MRCoSEAModeTrack(Modes, nextId, ...
                        lab, Xn, cfg);

                    %% Exploitation channel: manifold sampling of mature
                    %% modes replaces part of the SBX offspring (budget
                    %% neutral, 1:1; per-mode survival feedback adapts the
                    %% share, a mode whose model offspring keep dying
                    %% auto-degrades to pure SBX)
                    Offspring = [];
                    maniN     = 0;
                    modeTab   = zeros(0, 5);
                    modelCand = struct('id',{},'X',{});   % survival feedback
                    if cfg.manifoldEnabled
                        M    = Problem.M;
                        mTh  = max(1, min(M - 1, Problem.D - 1));
                        minS = max(cfg.minModeSize, 2*mTh + 2);
                        for k = 1 : numel(Modes)
                            if Modes(k).age < cfg.minModeAge || ...
                                    Modes(k).size < minS
                                continue;
                            end
                            mfrac = MRCoSEAManiFrac(fracA, Modes(k).id, cfg);
                            if mfrac <= 0
                                continue;
                            end
                            mem  = find(lab == k);
                            Xfit = Xn(mem, :);
                            fit  = MRCoSEAManifoldFit(Xfit, M, cfg);
                            if ~fit.valid
                                continue;
                            end
                            qModel = max(1, round(mfrac * Modes(k).size));
                            qModel = min(qModel, max(1, N - 2)); % keep GA alive
                            % v3.1 arm: absolute per-mode quota cap. The
                            % size-proportional share lets a degenerate
                            % giant mode claim 35% of the budget (Stage-0
                            % diagnosis: lastQ = 140 of N = 400), which
                            % amplifies the blob it feeds on.
                            qModel = min(qModel, ...
                                max(1, round(cfg.maniQCapFrac * N)));
                            cand = MRCoSEAManifoldSample(fit, qModel, ...
                                progress, cfg, Problem);
                            Offspring = [Offspring, Problem.Evaluation(cand)];
                            maniN = maniN + qModel;
                            modelCand(end+1) = struct( ...
                                'id',Modes(k).id,'X',{cand}); %#ok<AGROW>
                            modeTab(end+1, :) = [gen, Modes(k).id, ...
                                Modes(k).size, Modes(k).age, qModel]; %#ok<AGROW>
                        end
                    end

                    %% GA part: plain global binary tournament on the
                    %% mode-relative rank (mr on) or the global rank (mr
                    %% off, manifold-only ablation arm); unrestricted
                    %% mating restores cross-modal mixing (v2.0 pool
                    %% mating caused the G1 regression)
                    qGA  = N - maniN;
                    qGAe = 2 * floor(qGA / 2);    % OperatorGA halves parents
                    if qGAe > 0
                        if mrOn
                            rk = MRank;
                        else
                            rk = FrontNo;
                        end
                        mp = TournamentSelection(2, qGAe, rk, -Div);
                        Offspring = [Offspring, ...
                            OperatorGA(Problem, Population(mp))];
                    end

                    %% MR per-generation log
                    Algorithm.log.gen(end+1)    = gen;
                    Algorithm.log.nModes(end+1) = max([lab, 0]);
                    Algorithm.log.maniN(end+1)  = maniN;
                end

                %% Injection (reflection / random proposals / immigrants),
                %  every candidate is really evaluated (E5-safe, v0.8)
                if cfg.injectEnabled && injEnabled && gen >= startGen && ...
                        gen <= stopGen && mod(gen, kInj) == 0
                    m = min(mAdapt, maxFE - Problem.FE - N);
                    m = max(0, m);
                    if m > 0
                        if cfg.randomImmigrant
                            % A3 control (destructive): m evaluated random
                            % points replace m randomly chosen individuals
                            Ximm = unifrnd(repmat(Problem.lower,m,1), ...
                                           repmat(Problem.upper,m,1));
                            Imm  = Problem.Evaluation(Ximm);
                            kdel = min(m, length(Population));
                            delIdx = randperm(length(Population), kdel);
                            Population(delIdx) = [];
                            lastCand  = Ximm;
                            Offspring = [Offspring, Imm];
                        else
                            cand = CoSEAPropose(Population, Div, Problem, m, cfg, progress);
                            if ~isempty(cand)
                                lastCand  = cand;
                                Offspring = [Offspring, Problem.Evaluation(cand)];
                            end
                        end
                    end
                end

                %% Environmental selection: the mode-relative operator
                if mrOn
                    [Population, MRank, Div, info] = MRCoSEAEnvSelect( ...
                        [Population, Offspring], N, Problem, cfg, progress);
                    Algorithm.log.protN(end+1)  = info.protN;
                    Algorithm.log.maxModeSize(end+1) = info.maxModeSize;
                    Algorithm.log.floorSkip(end+1)   = info.floorSkip;
                    Algorithm.log.floorAdd(end+1)    = info.floorAdd;
                else
                    [Population, FrontNo, Div] = CoSEAEnvSelect( ...
                        [Population, Offspring], N, Problem, cfg, progress);
                end

                %% Per-mode survival feedback of the manifold channel:
                % model offspring pass selection unchanged, so exact
                % matching is valid (same argument as M3)
                if newOn && ~isempty(modelCand)
                    Xnow = Population.decs;         % SOLUTION.decs no index
                    for mc = 1 : numel(modelCand)
                        sM = CoSEASurvivalRate(modelCand(mc).X, Xnow);
                        fracA = MRCoSEAManiUpdate(fracA, ...
                            modelCand(mc).id, sM, cfg);
                    end
                end

                %% M3: survival-rate feedback and auto-degradation (v0.8)
                if ~isempty(lastCand)
                    Xnow = Population.decs;         % SOLUTION.decs takes no index
                    s = CoSEASurvivalRate(lastCand, Xnow);
                    Algorithm.log.inj  = [Algorithm.log.inj,  gen];
                    Algorithm.log.surv = [Algorithm.log.surv, s];
                    Algorithm.log.m    = [Algorithm.log.m,    size(lastCand,1)];
                    if cfg.adaptiveInj
                        if s > cfg.sHi
                            mAdapt = min(mMax, ceil(mAdapt * (1 + cfg.alpha)));
                        elseif s < cfg.sLo
                            mAdapt = max(mMin, floor(mAdapt / 2));
                        end
                        if s == 0
                            zeroCnt = zeroCnt + 1;
                            if zeroCnt >= cfg.zeroWin
                                injEnabled = false; % auto-degrade to backbone
                            end
                        else
                            zeroCnt = 0;
                        end
                    end
                    lastCand = [];
                end

                %% Mode table log (MR path rows, appended per generation)
                if newOn && ~isempty(modeTab)
                    Algorithm.log.modeTab = [Algorithm.log.modeTab; modeTab];
                end

                %% Checkpoints at every 20% of the budget (v0.8)
                if Problem.FE >= nextCkpt
                    Algorithm.ckpt{end+1} = {Problem.FE, Population.decs, Population.objs};
                    nextCkpt = nextCkpt + 0.2 * maxFE;
                end
            end
        end
    end
end

function f = MRCoSEAManiFrac(fracA, id, cfg)
%MRCoSEAManiFrac Effective model share of a mode (cap when untracked).
    if ~cfg.maniFeedback
        f = cfg.manifoldFrac;
        return;
    end
    hit = [fracA.id] == id;
    if any(hit)
        f = min(fracA(hit).f, cfg.manifoldFrac);
    else
        f = cfg.manifoldFrac;
    end
end

function fracA = MRCoSEAManiUpdate(fracA, id, s, cfg)
%MRCoSEAManiUpdate Survival feedback step of one mode's model share.
%   s < maniSLo -> halve (0 below 0.05 = per-mode auto-degradation to SBX);
%   s > maniSHi -> grow by 25% toward the configured share. Untouched
%   between the bands (hysteresis, same philosophy as M3).
    hit = [fracA.id] == id;
    if any(hit)
        ix = find(hit);
        f  = fracA(ix).f;
        if s < cfg.maniSLo
            f = f / 2;
            if f < 0.05
                f = 0;
            end
        elseif s > cfg.maniSHi
            f = min(1, f * 1.25);
        end
        fracA(ix).f = f;
    else
        f = cfg.manifoldFrac;                     % first observation
        if s < cfg.maniSLo
            f = f / 2;
        end
        fracA(end+1) = struct('id',id,'f',f); %#ok<AGROW>
    end
end
