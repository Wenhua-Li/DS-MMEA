function cfg = MRCoSEAConfig(userCfg)
%MRCoSEAConfig Default configuration of MRCoSEA (CoSEA v3.0) with override.
%   cfg = MRCoSEAConfig(userCfg) merges the user struct over the defaults.
%   All v0.8 (CoSEA) fields are inherited with identical default values, so
%   the configuration
%       struct('mrEnabled',false,'manifoldEnabled',false)
%   reproduces CoSEA v0.8 bit-for-bit (the off-path calls the v0.8
%   functions directly). MR-specific fields (v3.1 philosophy: ONE
%   selection operator, pure mode-relative survival; the v3.0 annealed
%   gate is deleted by evidence -- the Stage-1/2 races showed every
%   clustering and gate alternative is net-negative or dominated):
%     clusterMode/epsFrac                        - mode identification
%     gateType('none')                           - pure mode-relative
%     manifoldEnabled/...                        - per-mode local-PCA
%                                                manifold sampling
%   No workspace globals are used (parfor-safe).

    def = struct( ...
        ... % --- core mechanism switches (v0.8, unchanged) ---
        'injectEnabled',    true, ...
        'gateMode',         'off', ...
        'mFrac',            0.30, ...
        'kFrac',            0.05, ...
        'injectStartFrac',  0.10, ...
        'injectStopFrac',   1.00, ...
        ... % --- reflection proposal (v0.8, unchanged) ---
        'proposeMode',      'reflect', ...
        'epsFrac',          0.08, ...
        ... % --- A3 control (v0.8, unchanged) ---
        'randomImmigrant',  false, ...
        ... % --- environmental selection backbone (v0.8, unchanged) ---
        'envSel',           'speaEps', ...
        'noveltyFilter',    false, ...
        'injEpsFrac',       0.15, ...
        ... % --- M3 survival feedback (v0.8, unchanged) ---
        'adaptiveInj',      true, ...
        'alpha',            0.25, ...
        'sHi',              0.20, ...
        'sLo',              0.02, ...
        'mMinFrac',         0.05, ...
        'zeroWin',          3, ...
        ... % --- M4 reflection quality (v0.8, unchanged) ---
        'noiseDecay',       true, ...
        'sigma0',           0.02, ...
        'sigmaEnd',         0.005, ...
        'srcWeight',        true, ...
        'minClusterQuota',  2, ...
        ... % --- M8 retention (v0.8, unchanged; used only on the
        ... %     mrEnabled=false paths as the v0.8 selection) ---
        'epsRetain0',       0.30, ...
        'epsRetainMin',     0.02, ...
        'epsRetainCapFrac', 0.10, ...
        'mateRestrict',     false, ...
        ... % --- MR v3.0: mode identification (v2.0 semantics) ---
        'clusterMode',      'eps', ...    % 'eps' radius | 'knn' graph
        'knnK',             [], ...       % k of knn graph ([] = round(sqrt(n)))
        'matchThrFrac',     1.0, ...      % match radius = frac*epsFrac*sqrt(d)
        ... % --- MR v3.0: mode-relative selection ---
        'mrEnabled',        true, ...     % the one selection operator
        'epsAnneal0',       0.60, ...     % gate epsilon at progress 0
        'epsFloor',         0.25, ...     % gate epsilon floor (late budget)
        'modeFloor',        0, ...        % ungated per-mode front-1 floor (OFF by
        ...                              % default: the P3.2 experiments showed
        ...                              % it is either inert or chaos-dominated
        ...                              % -- the merged population always holds
        ...                              % one ~2N blob, so the size-guarded
        ...                              % floor adds ~0 members; kept as an
        ...                              % ablation arm only)
        'modeFloorMaxFrac', 0.5, ...      % floor applies only to modes with
        ...                              % size <= frac*N (developing modes;
        ...                              % the degenerate single-blob regime
        ...                              % is excluded)
        'annealShape',      'commea', ... % 'commea' (-log2(1.5t): loose while
        ...                              % modes develop, tightens late; won
        ...                              % the P1 smoke on G2 developing
        ...                              % modes 5x over 'linear')
        ...                              % 'linear' kept as the ablation arm
        ... % --- MR v3.0: per-mode manifold learning (exploitation) ---
        'manifoldEnabled',  true, ...     % model share cap (mature modes)
        'manifoldFrac',     0.35, ...
        'minModeSize',      8, ...        % min members before modelling
        'minModeAge',       3, ...        % min tracked generations
        'varThresh',        0.90, ...     % cumulative-variance latent cap
        'residMax',         0.20, ...     % residual-variance validity guard
        'maniNoise0',       1.0, ...      % residual noise multiplier (decays)
        'maniNoiseEnd',     0.25, ...
        'maniFeedback',     true, ...     % per-mode survival-rate feedback of
        'maniSLo',          0.02, ...     % model offspring (auto-degrade to
        'maniSHi',          0.20, ...     % SBX per mode when it stops paying)
        ... % --- MR v3.1 smoke arms (defaults reproduce v3.0 exactly) ---
        ... % radMult: 'quantile'/'pca' cluster radius = radMult * the
        ... % nnQuantile percentile of the 1-NN distance distribution
        ... % (data-adaptive, replaces the epsFrac*sqrt(d) absolute length
        ... % that chains high-D populations into blobs).
        'radMult',          2.0, ...
        ... % nnQuantile: which percentile of the 1-NN distances sets the
        ... % radius. 50 (median) = Stage-1 'qrad' arm, shatters dense low-D
        ... % PS curves at sampling gaps; 90 (Q90, 'q90' arm) is gap-robust:
        ... % the radius then tracks the spacing TAIL, so curves stay whole
        ... % in low-D while high-D radii shrink ~25x (offline check: true
        ... % mode counts recovered on Polygon/IDMP, K=1-4 kept on MMF).
        'nnQuantile',       50, ...
        ... % maniQCapFrac: per-mode per-generation model quota cap as a
        ... % fraction of N; 1.0 = no cap (v3.0). Guards the rich-get-
        ... % richer blob amplification found in the Stage-0 diagnosis.
        'maniQCapFrac',     1.0, ...
        ... % gateType: v3.1 DEFAULT 'none' = pure mode-relative selection
        ... % (no global pruning at all). The Stage-1/2 races (2026-09-28/29,
        ... % 12 arms) showed every gate variant is dominated by deletion:
        ... % the per-member global gate fixes no rank on average (P4:
        ... % p = 0.89) but executes the weak-mode colonies of G2-M4/G5
        ... % (17x/11x disasters). 'global' keeps the v3.0 annealed
        ... % per-member gate as an ablation arm; 'modeNet'/'dualRef' are
        ... % the rejected race variants (kept for the record).
        'gateType',         'none', ...
        ... % robustNorm: gate normalization base = the global front-1
        ... % range (true) instead of the merged min/max (false = v3.0),
        ... % so objective-space outliers cannot stretch the eps scale.
        'robustNorm',       false);

    cfg = def;
    if isstruct(userCfg)
        f = fieldnames(userCfg);
        for i = 1 : numel(f)
            cfg.(f{i}) = userCfg.(f{i});
        end
    end
end
