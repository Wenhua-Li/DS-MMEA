function cfg = CoSEAConfig(userCfg)
%CoSEAConfig Default configuration of CoSEA with user override (M1).
%   cfg = CoSEAConfig(userCfg) merges the user struct over the defaults.
%   Call from within CoSEA.main as CoSEAConfig(Algorithm.ParameterSet(struct()))
%   (ParameterSet is protected, so the call must live in a class method).
%   No workspace globals are used (parfor workers receive cfg explicitly).

    def = struct( ...
        ... % --- core mechanism switches ---
        'injectEnabled',    true, ...    % false = pure backbone (degenerate test)
        'gateMode',         'off', ...   % cluster-growth gate (kept off, falsified)
        'mFrac',            0.30, ...    % injected candidates per injection = mFrac*N
        'kFrac',            0.05, ...    % injection interval = kFrac*maxgen
        'injectStartFrac',  0.10, ...    % no injection in the first 10% generations
        'injectStopFrac',   1.00, ...    % no injection after this budget fraction (A5)
        ... % --- reflection proposal ---
        'proposeMode',      'reflect', ... % 'reflect' (main) | 'random' (A3 control)
        'epsFrac',          0.08, ...    % cluster radius = epsFrac*sqrt(d)
        ... % --- A3 control: random immigrants (destructive restart) ---
        'randomImmigrant',  false, ...   % replace m individuals by m random points
        ... % --- environmental selection ---
        'envSel',           'speaEps', ... % 'crowding'(v0.6) | 'spea'(v0.7) | 'speaEps'(v0.8/M8)
        ... % --- novelty filter (optional, off) ---
        'noveltyFilter',    false, ...
        'injEpsFrac',       0.15, ...
        ... % --- M3: survival-rate feedback, adaptive injection, auto-degrade ---
        'adaptiveInj',      true, ...
        'alpha',            0.25, ...    % multiplicative growth step of mAdapt
        'sHi',              0.20, ...    % survival rate above this: increase m
        'sLo',              0.02, ...    % survival rate below this: halve m
        'mMinFrac',         0.05, ...    % lower bound of injection size (of N)
        'zeroWin',          3, ...       % consecutive zero-survival injections -> off
        ... % --- M4: reflection quality scheduling ---
        'noiseDecay',       true, ...    % sigma: sigma0 -> sigmaEnd linearly over budget
        'sigma0',           0.02, ...
        'sigmaEnd',         0.005, ...
        'srcWeight',        true, ...    % sample reflection sources weighted by isolation
        'minClusterQuota',  2, ...       % per-cluster minimum quota (small-mode guard)
        ... % --- M8: local-PS retention (epsilon-dominance) + mating restriction ---
        'epsRetain0',       0.30, ...    % initial multiplicative epsilon
        'epsRetainMin',     0.02, ...    % lower bound of epsilon (late stage)
        'epsRetainCapFrac', 0.10, ...    % retained set size cap (fraction of N)
        'mateRestrict',     false);      % decision-space nearest-neighbor mating

    cfg = def;
    if isstruct(userCfg)
        f = fieldnames(userCfg);
        for i = 1 : numel(f)
            cfg.(f{i}) = userCfg.(f{i});
        end
    end
end
