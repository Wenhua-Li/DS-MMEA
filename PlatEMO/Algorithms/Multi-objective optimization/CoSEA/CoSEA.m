classdef CoSEA < ALGORITHM
% <2026> <multi> <real> <multimodal>
% CoSEA - coverage-aware evolutionary algorithm for multimodal multiobjective
%   optimization, rebuilt 2026-09 on the class-based PlatEMO (v0.8):
%     backbone  - NSGA-II non-dominated sorting + decision-space diversity
%                 maintenance (CoSEAEnvSelect: crowding / spea / speaEps)
%     discovery - manifold-aligned symmetric reflection (CoSEAReflectPropose)
%   Injected candidates are really evaluated and counted in the budget
%   (Problem.Evaluation), no free evaluations.
%   M3: survival-rate feedback adapts the injection size and auto-degrades
%   to the pure backbone when injected candidates never survive (e.g. on
%   problems without central symmetry).
%   M8: local-PS retention via multiplicative epsilon-dominance protects
%   newborn/rare modes from convergence pressure.
%
%   Configuration: CoSEA('parameter',{struct('envSel','spea',...)}) overrides
%   any subset of CoSEAConfig defaults (no workspace globals, parfor-safe).

%------------------------------- Copyright --------------------------------
% Reconstruction of the CoSEA algorithm (thesis project, 2026TEVC S2).
%--------------------------------------------------------------------------

    properties
        ckpt = {};  % checkpoints: {FE, PopDec, PopObj} at every 20% budget
        log  = struct('inj',[],'surv',[],'m',[]); % injection/survival log
    end
    methods
        function main(Algorithm, Problem)
            cfg = CoSEAConfig(Algorithm.ParameterSet(struct()));
            N     = Problem.N;
            maxFE = Problem.maxFE;

            %% Scheduling constants (in generations)
            maxgen   = ceil((maxFE - N) / N);
            kInj     = max(2, round(cfg.kFrac * maxgen));
            startGen = ceil(cfg.injectStartFrac * maxgen);
            % A5 (injection timing): stopGen = Inf keeps the v0.8 behaviour
            % bit-identical, since gen never exceeds maxgen.
            if isfield(cfg,'injectStopFrac') && cfg.injectStopFrac < 1
                stopGen = floor(cfg.injectStopFrac * maxgen);
            else
                stopGen = Inf;
            end

            %% M3 state
            mMax       = round(cfg.mFrac * N);
            mMin       = max(1, round(cfg.mMinFrac * N));
            mAdapt     = mMax;
            injEnabled = true;
            zeroCnt    = 0;
            lastCand   = [];

            %% Generate random population
            Population = Problem.Initialization();
            [Population, FrontNo, Div] = CoSEAEnvSelect(Population, N, Problem, cfg, 0);
            gen = 0;
            nextCkpt = 0.2 * maxFE;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                progress = Problem.FE / maxFE;
                MatingPool = CoSEAMateSelect(Population, FrontNo, Div, Problem, cfg);
                Offspring  = OperatorGA(Problem, Population(MatingPool));

                %% Injection (reflection / random proposals / immigrants),
                %  every candidate is really evaluated (E5-safe)
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
                            k    = min(m, length(Population));
                            Population(randperm(length(Population),k)) = [];
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

                %% Environmental selection (backbone)
                [Population, FrontNo, Div] = CoSEAEnvSelect([Population, Offspring], N, Problem, cfg, progress);

                %% M3: survival-rate feedback and auto-degradation
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

                %% Checkpoints at every 20% of the budget
                if Problem.FE >= nextCkpt
                    Algorithm.ckpt{end+1} = {Problem.FE, Population.decs, Population.objs};
                    nextCkpt = nextCkpt + 0.2 * maxFE;
                end
            end
        end
    end
end
