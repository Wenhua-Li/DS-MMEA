classdef MSCoSEA < ALGORITHM
% <2026> <multi> <real> <multimodal>
% MSCoSEA - Multi-Scale CoSEA (CoSEA v4.0 candidate): one selection
% principle (mode-relative survival, MRCoSEA v3.1) running on K parallel
% clustering-scale hypotheses.
%   The 12-arm race of 2026-09-28/29 closed the single-population design
%   space: curve-like modes (G1/G3) need a loose clustering radius, compact
%   high-dimensional modes (G4/G6) need a tight one, and no single rule can
%   be right on both (the correct scale is a property of the problem
%   geometry that on-track information cannot infer reliably). MSCoSEA does
%   not infer the scale -- it runs two structural hypotheses in parallel
%   and lets survival arbitrate between them:
%     family W (wide) - the v3.1 configuration verbatim ('eps', 0.08*sqrt(d)
%                radius): correct on curves and sparse developing colonies.
%     family T (tight) - the same selection operator under the 'epsq'
%                clustering (Stage-2 calibration 2xQ90 1-NN cap): separates
%                equivalent modes where the wide radius fuses them into one
%                blob (the G2-M4 bistability and G4/G6 fusion failures).
%                Its low-D fragmentation is harmless by construction: the
%                fragments' within-mode front-1 members are still good
%                solutions and the wide family carries the development.
%   The two families share ONE offspring pool per generation (CoMMEA's
%   coevolution pattern): GA offspring (global binary tournament over the
%   merged 2N parents), manifold samples of every family's mature modes,
%   and the reflection injections are evaluated ONCE and offered to both
%   environmental selections. Budget accounting is unchanged from v3.1:
%   per generation the total evaluation is ~N + injections, family
%   duplication costs no function evaluations (the initial population is
%   copied, not re-evaluated).
%   The output population is the joint SPEA2 truncation of all families to
%   N (display/checkpoint/metric only, never fed back into the search).
%   nPop = 1 reduces the loop onto the exact MRCoSEA v3.1 call sequence
%   (regression-gated bit-identity by exp/ms/test_ms_equiv.m).
%   New active fields: 2 (nPop, clusterModeB) -- see MSCoSEAConfig.

%------------------------------- Copyright --------------------------------
% Thesis project 2026TEVC S2 (CoSEA v4.0 multi-scale revision, 2026-09).
%--------------------------------------------------------------------------

    properties
        ckpt = {};  % checkpoints: {FE, PopDec, PopObj} at every 20% budget
        log  = struct('inj',[],'surv',[],'m',[], ...      % shared injection log
            'gen',[],'nModesK',[],'maniN',[], ... % per-generation, families
            'protNK',[],'maxModeSizeK',[], ...    % per-family selection info
            'modeTab',zeros(0,6)); % [gen fam id size age q]
        snap = {};  % opt-in manifold-channel snapshots (cfg.maniSnap):
                    % one struct per due generation, fields tag/gen/fe/
                    % progress/X/obj/lab/modesF/fracAF/cand (record-only,
                    % see localSnapRec)
    end
    methods
        function main(Algorithm, Problem)
            cfg = MSCoSEAConfig(Algorithm.ParameterSet(struct()));
            if cfg.segMani && cfg.quadMani
                error('MSCoSEA:maniArms', ...
                    'segMani and quadMani are mutually exclusive race arms');
            end
            N     = Problem.N;
            maxFE = Problem.maxFE;
            K     = max(1, round(cfg.nPop));

            %% Per-family configs: family 1 = v3.1 verbatim (wide); family 2
            %% = the tight hypothesis; family 3+ = the C overrides (inherit
            %% family-2 when empty; the K = 3 three-scale arm)
            cfgF = cell(1, K);
            for k = 1 : K
                ck = cfg;
                if k >= 2
                    ck.clusterMode = cfg.clusterModeB;
                    if ~isempty(cfg.epsFracB)
                        ck.epsFrac = cfg.epsFracB;
                    end
                    if ~isempty(cfg.radMultB)
                        ck.radMult = cfg.radMultB;
                    end
                    if ~isempty(cfg.nnQuantileB)
                        ck.nnQuantile = cfg.nnQuantileB;
                    end
                end
                if k >= 3
                    if ~isempty(cfg.clusterModeC)
                        ck.clusterMode = cfg.clusterModeC;
                    end
                    if ~isempty(cfg.epsFracC)
                        ck.epsFrac = cfg.epsFracC;
                    end
                    if ~isempty(cfg.radMultC)
                        ck.radMult = cfg.radMultC;
                    end
                    if ~isempty(cfg.nnQuantileC)
                        ck.nnQuantile = cfg.nnQuantileC;
                    end
                end
                cfgF{k} = ck;
            end

            %% Scheduling constants (in generations, v0.8 semantics)
            maxgen   = ceil((maxFE - N) / N);
            kInj     = max(2, round(cfg.kFrac * maxgen));
            startGen = ceil(cfg.injectStartFrac * maxgen);
            if isfield(cfg,'injectStopFrac') && cfg.injectStopFrac < 1
                stopGen = floor(cfg.injectStopFrac * maxgen);
            else
                stopGen = Inf;
            end

            %% M3 state (v0.8): ONE shared discovery channel for all families
            mMax       = round(cfg.mFrac * N);
            mMin       = max(1, round(cfg.mMinFrac * N));
            mAdapt     = mMax;
            injEnabled = true;
            zeroCnt    = 0;
            lastCand   = [];

            %% Mode state per family (independent id spaces): persistent
            %% lists for age tracking and the per-mode model-share feedback
            ModesF  = cell(1, K);
            fracAF  = cell(1, K);
            nextIdF = zeros(1, K);
            for k = 1 : K
                ModesF{k} = struct('id',{},'cent',{},'size',{},'age',{}); %#ok<STRUCT>
                fracAF{k} = struct('id',{},'f',{});          % adaptive shares
            end

            %% Generate random population; copy it to every family WITHOUT
            %% re-evaluation (FE starts at N exactly like v3.1)
            Population = Problem.Initialization();
            Pop   = cell(1, K);
            MRank = cell(1, K);
            Div   = cell(1, K);
            for k = 1 : K
                Pop{k} = Population;
            end
            for k = 1 : K
                [Pop{k}, MRank{k}, Div{k}] = MRCoSEAEnvSelect( ...
                    Pop{k}, N, Problem, cfgF{k}, 0);
            end
            OutPop   = MSCoSEAFuse(Pop, cfgF, N, Problem);
            gen      = 0;
            nextCkpt = 0.2 * maxFE;

            %% Optimization
            while Algorithm.NotTerminated(OutPop)
                gen = gen + 1;
                progress = Problem.FE / maxFE;

                %% Track modes on the parents, per family, at its own scale
                lab  = cell(1, K);
                XnF  = cell(1, K);
                for k = 1 : K
                    X       = Pop{k}.decs;        % E2: no indexed .decs
                    lab{k}  = MRCoSEAModeCluster(X, Problem, cfgF{k});
                    XnF{k}  = (X - Problem.lower) ./ ...
                              max(Problem.upper - Problem.lower, eps);
                    [ModesF{k}, nextIdF(k)] = MRCoSEAModeTrack( ...
                        ModesF{k}, nextIdF(k), lab{k}, XnF{k}, cfgF{k});
                end

                %% Manifold snapshot (opt-in logging, default off; added
                %% 2026-09-30 for the mechanism-visualization report).
                %% Records the exact family state this generation's
                %% manifold channel sees; the actual model offspring are
                %% attached below the channel. Pure recording -- no RNG
                %% calls, no search interaction. With maniSnap = false
                %% the branch never executes (T1' bit-identity intact).
                snapIdx = 0;
                if cfg.maniSnap
                    snapIdx = localSnapRec(Algorithm, gen, Problem.FE, ...
                        maxFE, N, Pop, lab, ModesF, fracAF);
                end

                %% Exploitation channel: manifold sampling of every
                %% family's mature modes (budget-neutral 1:1; per-mode
                %% survival feedback adapts the share). On the tight family
                %% fragmented modes stay immature and are skipped for free.
                Offspring = [];
                maniN     = 0;
                modeTab   = zeros(0, 6);
                modelCand = struct('fam',{},'id',{},'X',{});
                if cfg.manifoldEnabled
                    M    = Problem.M;
                    mTh  = max(1, min(M - 1, Problem.D - 1));
                    minS = max(cfg.minModeSize, 2*mTh + 2);
                    poolFull = false;
                    for k = 1 : K
                        if poolFull
                            break;
                        end
                        for j = 1 : numel(ModesF{k})
                            md = ModesF{k}(j);
                            if md.age < cfg.minModeAge || ...
                                    md.size < minS
                                continue;
                            end
                            mfrac = MSCoSEAManiFrac(fracAF{k}, md.id, cfg);
                            if mfrac <= 0
                                continue;
                            end
                            mem  = find(lab{k} == j);
                            if cfg.segMani || cfg.quadMani
                                %% Manifold-upgrade race arms (2026-10-01,
                                %% both default false; the else branch below
                                %% is the v4.1 incumbent verbatim, so the
                                %% default path and T1' bit-identity are
                                %% untouched). Quota caps identical in all
                                %% arms; only the model class differs.
                                if cfg.segMani
                                    segs = MRCoSEAManifoldFitSeg( ...
                                        XnF{k}(mem, :), M, minS, cfg);
                                    if isempty(segs)
                                        continue;
                                    end
                                    qTot = max(1, round(mfrac * md.size));
                                    qTot = min(qTot, max(1, N - 2));
                                    qTot = min(qTot, ...
                                        max(1, round(cfg.maniQCapFrac * N)));
                                    qTot = min(qTot, max(0, N - 2 - maniN));
                                    if qTot < 1
                                        poolFull = true;
                                        break;
                                    end
                                    cand  = [];
                                    cumSz = cumsum([segs.size]);
                                    prev  = 0;
                                    for s = 1 : numel(segs)
                                        qs = round(cumSz(s) / md.size * qTot) - ...
                                             round(prev / md.size * qTot);
                                        prev = cumSz(s);
                                        if qs < 1
                                            continue;
                                        end
                                        csamp = MRCoSEAManifoldSample( ...
                                            segs(s).fit, qs, ...
                                            progress, cfg, Problem);
                                        cand = [cand; csamp]; %#ok<AGROW>
                                    end
                                    if isempty(cand)
                                        continue;
                                    end
                                    Offspring = [Offspring, ...
                                        Problem.Evaluation(cand)];
                                    maniN = maniN + size(cand, 1);
                                    modelCand(end+1) = struct( ...
                                        'fam',k,'id',md.id,'X',{cand}); %#ok<AGROW>
                                    modeTab(end+1, :) = [gen, k, md.id, ...
                                        md.size, md.age, size(cand,1)]; %#ok<AGROW>
                                else
                                    fq = MRCoSEAManifoldFitQuad( ...
                                        XnF{k}(mem, :), M, cfg);
                                    if ~fq.valid
                                        continue;
                                    end
                                    qModel = max(1, round(mfrac * md.size));
                                    qModel = min(qModel, max(1, N - 2));
                                    qModel = min(qModel, ...
                                        max(1, round(cfg.maniQCapFrac * N)));
                                    qModel = min(qModel, ...
                                        max(0, N - 2 - maniN));
                                    if qModel < 1
                                        poolFull = true;
                                        break;
                                    end
                                    cand = MRCoSEAManifoldSampleQuad(fq, ...
                                        qModel, progress, cfg, Problem);
                                    Offspring = [Offspring, ...
                                        Problem.Evaluation(cand)];
                                    maniN = maniN + qModel;
                                    modelCand(end+1) = struct( ...
                                        'fam',k,'id',md.id,'X',{cand}); %#ok<AGROW>
                                    modeTab(end+1, :) = [gen, k, md.id, ...
                                        md.size, md.age, qModel]; %#ok<AGROW>
                                end
                            else
                                fit  = MRCoSEAManifoldFit(XnF{k}(mem, :), M, cfg);
                                if ~fit.valid
                                    continue;
                                end
                                qModel = max(1, round(mfrac * md.size));
                                qModel = min(qModel, max(1, N - 2)); % keep GA alive
                                qModel = min(qModel, ...
                                    max(1, round(cfg.maniQCapFrac * N)));
                                % Cross-family fuse: once the model share would
                                % starve GA, stop. Never triggers for K = 1
                                % (sum of per-mode shares <= 0.35N there).
                                qModel = min(qModel, max(0, N - 2 - maniN));
                                if qModel < 1
                                    poolFull = true;
                                    break;
                                end
                                cand = MRCoSEAManifoldSample(fit, qModel, ...
                                    progress, cfg, Problem);
                                Offspring = [Offspring, Problem.Evaluation(cand)];
                                maniN = maniN + qModel;
                                modelCand(end+1) = struct( ...
                                    'fam',k,'id',md.id,'X',{cand}); %#ok<AGROW>
                                modeTab(end+1, :) = [gen, k, md.id, ...
                                    md.size, md.age, qModel]; %#ok<AGROW>
                            end
                        end
                    end
                end

                if snapIdx > 0
                    S = Algorithm.snap{snapIdx};
                    S.cand = modelCand;
                    Algorithm.snap{snapIdx} = S;
                end

                %% GA part: plain global binary tournament on the MERGED
                %% parents (all families, cross-family parities). A solution
                %% that is mode-front-1 under EITHER scale hypothesis has
                %% equal standing -- no quotas, no family-restricted mating.
                qGA  = N - maniN;
                qGAe = 2 * floor(qGA / 2);    % OperatorGA halves parents
                if qGAe > 0
                    allP  = [Pop{:}];
                    rk    = [MRank{:}];
                    dv    = [Div{:}];
                    mp    = TournamentSelection(2, qGAe, rk, -dv);
                    Offspring = [Offspring, OperatorGA(Problem, allP(mp))];
                end

                %% Per-generation log
                Algorithm.log.gen(end+1)      = gen;
                Algorithm.log.nModesK(end+1,1:K) = ...
                    cellfun(@(l) max([l, 0]), lab);
                Algorithm.log.maniN(end+1)    = maniN;

                %% Injection (reflection / random proposals / immigrants):
                %% one channel, proposed from the merged families, every
                %% candidate really evaluated (E5-safe, v0.8)
                if cfg.injectEnabled && injEnabled && gen >= startGen && ...
                        gen <= stopGen && mod(gen, kInj) == 0
                    m = min(mAdapt, maxFE - Problem.FE - N);
                    m = max(0, m);
                    if m > 0
                        if cfg.randomImmigrant
                            % A3 control (destructive): m evaluated random
                            % points replace m randomly chosen members of
                            % family 1 (the wide family). K = 1 reproduces
                            % the v0.8/MRCoSEA semantics exactly.
                            Ximm = unifrnd(repmat(Problem.lower,m,1), ...
                                           repmat(Problem.upper,m,1));
                            Imm  = Problem.Evaluation(Ximm);
                            kdel = min(m, length(Pop{1}));
                            delIdx = randperm(length(Pop{1}), kdel);
                            Pop{1}(delIdx) = [];
                            lastCand  = Ximm;
                            Offspring = [Offspring, Imm];
                        else
                            allP  = [Pop{:}];
                            allDv = [Div{:}];
                            cand = CoSEAPropose(allP, allDv, ...
                                Problem, m, cfg, progress);
                            if ~isempty(cand)
                                lastCand  = cand;
                                Offspring = [Offspring, Problem.Evaluation(cand)];
                            end
                        end
                    end
                end

                %% Environmental selection: the SAME mode-relative operator,
                %% once per family, each on [own family + (shared or sliced)
                %% offspring]. shareOffspring = false is the isolated-islands
                %% ablation: the pool is generated once as usual but family k
                %% only receives the offspring positions interleaved to it,
                %% so one family's search products no longer flow to the
                %% other (mating/injection generation stays shared).
                infoF = cell(1, K);
                for k = 1 : K
                    if cfg.shareOffspring || isempty(Offspring)
                        kOff = Offspring;
                    else
                        kOff = Offspring(mod(0:length(Offspring)-1, K) == k-1);
                    end
                    [Pop{k}, MRank{k}, Div{k}, infoF{k}] = ...
                        MRCoSEAEnvSelect([Pop{k}, kOff], ...
                        N, Problem, cfgF{k}, progress);
                end
                Algorithm.log.protNK(end+1, 1:K) = ...
                    cellfun(@(s) s.protN, infoF);
                Algorithm.log.maxModeSizeK(end+1, 1:K) = ...
                    cellfun(@(s) s.maxModeSize, infoF);

                %% Per-mode survival feedback of the manifold channel:
                %% survival in ANY family counts (model offspring pass
                %% selection unchanged, so exact matching is valid)
                if ~isempty(modelCand)
                    allP  = [Pop{:}];
                    Xnow  = allP.decs;          % SOLUTION.decs no index
                    for mc = 1 : numel(modelCand)
                        sM = CoSEASurvivalRate(modelCand(mc).X, Xnow);
                        f  = modelCand(mc).fam;
                        fracAF{f} = MSCoSEAManiUpdate(fracAF{f}, ...
                            modelCand(mc).id, sM, cfg);
                    end
                end

                %% M3: survival-rate feedback and auto-degradation (v0.8);
                %% survival in ANY family counts
                if ~isempty(lastCand)
                    allP  = [Pop{:}];
                    Xnow  = allP.decs;         % SOLUTION.decs takes no index
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

                %% Mode table log (per-family rows, appended per generation)
                if ~isempty(modeTab)
                    Algorithm.log.modeTab = [Algorithm.log.modeTab; modeTab];
                end

                %% Output fusion: report-layer mode-front filter + tiered
                %% maximin truncation of all families to N (display/
                %% checkpoint population only, never fed back). K = 1
                %% returns the family untouched (bit-identity).
                OutPop = MSCoSEAFuse(Pop, cfgF, N, Problem);

                %% Checkpoints at every 20% of the budget (v0.8)
                if Problem.FE >= nextCkpt
                    Algorithm.ckpt{end+1} = {Problem.FE, OutPop.decs, OutPop.objs};
                    nextCkpt = nextCkpt + 0.2 * maxFE;
                end
            end
        end
    end
end

function idx = localSnapRec(Algorithm, gen, fe, maxFE, N, Pop, lab, ...
    ModesF, fracAF)
%localSnapRec Append one manifold-channel snapshot when due (record-only).
%   Schedule: 'init' at gen 1, then the first generation crossing 10% /
%   50% / 85% of the budget, plus 'final' at the last generation
%   (fe + N >= maxFE: the last generation whose evaluation still fits).
%   At most one snapshot per generation (priority final > late > mid >
%   early > init). Returns the snap index, 0 when no snapshot was due.

    K = numel(Pop);
    if isempty(Algorithm.snap)
        used = {};
    else
        used = cellfun(@(s) s.tag, Algorithm.snap, 'UniformOutput', false);
    end
    tag = '';
    if gen == 1
        tag = 'init';
    elseif fe + N >= maxFE
        tag = 'final';
    elseif fe / maxFE >= 0.85
        tag = 'late';
    elseif fe / maxFE >= 0.50
        tag = 'mid';
    elseif fe / maxFE >= 0.10
        tag = 'early';
    end
    idx = 0;
    if isempty(tag) || any(strcmp(used, tag))
        return;
    end
    S = struct();
    S.tag      = tag;
    S.gen      = gen;
    S.fe       = fe;
    S.progress = fe / maxFE;
    S.X   = cell(1, K);                      % SOLUTION.decs: whole-array
    S.obj = cell(1, K);                      % access only (env pit E2)
    for k = 1 : K
        S.X{k}   = Pop{k}.decs;
        S.obj{k} = Pop{k}.objs;
    end
    S.lab    = lab;
    S.modesF = ModesF;
    S.fracAF = fracAF;
    S.cand   = struct('fam', {}, 'id', {}, 'X', {});
    Algorithm.snap{end+1} = S;
    idx = numel(Algorithm.snap);
end

function OutPop = MSCoSEAFuse(Pop, cfgF, N, Problem)
%MSCoSEAFuse Joint output population of all families, truncated to N.
%   K = 1 (or a total <= N) returns the concatenation untouched (this fast
%   path precedes every filtering statement, keeping the nPop = 1
%   bit-identity gate intact). Otherwise the merged pool is clustered once
%   at the WIDE family scale and the delivered set is assembled by
%   MSCoSEAReport (2026-09-30 design revision): eligibility is
%   mode-relative where a mode carries information (clusters of at least
%   3 members: within-mode front-1) and pool-relative where it does not
%   (unclustered members and trivial 1-2 member clusters: front of the
%   merged pool); then greedy maximin over the eligible pool, refilled in
%   ascending pool-front order with accumulated nearest-distance when
%   fewer than N are eligible. The mode-relative branch is never global:
%   on benchmarks whose truth includes local Pareto fronts, globally
%   dominated members are load-bearing (71-problem deletion study,
%   2026-09-30). A dual-scale union variant (prune only if dominated at
%   both the wide and the tight scale) was panel-tested the same day and
%   rejected: it eliminates the rare single-run coverage tails but keeps
%   most of the within-mode-dominated residue the filter exists to remove
%   (panel record: Data_PF2 vs Data_PF3, exp/rq3). Zero new parameters;
%   the filter never feeds back into the search.

    allP = [Pop{:}];
    if numel(Pop) == 1 || length(allP) <= N
        OutPop = allP;
        return;
    end
    Xn   = (allP.decs - Problem.lower) ./ ...
           max(Problem.upper - Problem.lower, eps);
    labR = { MRCoSEAModeCluster(allP.decs, Problem, cfgF{1}) };
    sel  = MSCoSEAReport(Xn, allP.objs, allP.cons, labR, N);
    OutPop = allP(sel);
end

function f = MSCoSEAManiFrac(fracA, id, cfg)
%MSCoSEAManiFrac Effective model share of a mode (cap when untracked).
%   Identical semantics to the MRCoSEA local function (self-contained
%   copy under the unique-name rule).

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

function fracA = MSCoSEAManiUpdate(fracA, id, s, cfg)
%MSCoSEAManiUpdate Survival feedback step of one mode's model share.
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
