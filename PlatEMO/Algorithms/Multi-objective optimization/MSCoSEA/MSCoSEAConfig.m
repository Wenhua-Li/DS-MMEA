function cfg = MSCoSEAConfig(userCfg)
%MSCoSEAConfig Default configuration of MSCoSEA (CoSEA v4.0 candidate).
%   cfg = MSCoSEAConfig(userCfg) starts from the full MRCoSEA (v3.1) default
%   configuration (family 1 runs it unchanged: the wide family IS v3.1) and
%   adds the multi-scale fields:
%     nPop       - number of parallel scale families (default 2). nPop = 1
%                  collapses the main loop onto the single-family MRCoSEA
%                  v3.1 call sequence (regression gate T1' checks
%                  bit-identity against MRCoSEA on the same seed).
%     clusterModeB / epsFracB / radMultB / nnQuantileB - clustering of the
%                  tight families (2..K). clusterModeB defaults to 'epsq'
%                  (radius = min(epsFrac*sqrt(d), radMultB * the
%                  nnQuantileB percentile of 1-NN distances)), the Stage-2
%                  race arm whose high-D repair face (Polygon M4_D8 0.37x,
%                  M6_D8 0.40x, IDMPT1_NP4 0.03x) is inherited by the tight
%                  family while its low-D fragmentation harm is neutralized
%                  by the wide family. nnQuantileB defaults to 90 (the
%                  Stage-2 calibration 2xQ90; the v3.1 default 50 belongs
%                  to the rejected 'qrad' arm). The other B-fields default
%                  to [] = inherit the family-1 value (inert overrides for
%                  sensitivity arms only).
%   Active new fields: 2 (nPop, clusterModeB). No workspace globals
%   (parfor-safe).

    def = MRCoSEAConfig(struct());
    def.mrEnabled  = true;  % forced semantics: MSCoSEA has no global-selection
    ...               % path; that ablation arm lives in MRCoSEA (_global)
    def.nPop        = 2;
    def.clusterModeB = 'epsq';
    def.epsFracB    = [];
    def.radMultB    = [];
    def.nnQuantileB = 90;
    ... % shareOffspring: false = the isolated-islands ablation arm. The
    ... % offspring pool is still generated ONCE (merged mating, shared
    ... % reflection channel, budget unchanged) but family k only sees the
    ... % offspring slice interleaved to it (odd/even positions for K = 2),
    ... % so one family's search products no longer flow to the other.
    def.shareOffspring = true;
    ... % Family-3+ overrides (K = 3 three-scale arm; [] = inherit family-2).
    def.clusterModeC = [];
    def.epsFracC     = [];
    def.radMultC     = [];
    def.nnQuantileC  = [];
    ... % maniSnap: opt-in recording of manifold-channel snapshots for the
    ... % mechanism-visualization report (exp/ms/snap_mani_run.m). Pure
    ... % logging (populations, mode labels, fitted models, model
    ... % offspring); never affects the search. Default false.
    def.maniSnap = false;
    ... % Manifold-upgrade race arms (2026-10-01, dev plan S9; mutually
    ... % exclusive, both default false = v4.1 incumbent bit-identical):
    ... %   segMani  - segmented chords: recursive median split of a mode
    ... %              into pieces whose linear fit explains >= varThresh
    ... %              (MRCoSEAManifoldFitSeg; split bar = 1 - varThresh,
    ... %              floor = the mode-size gate, piece cap 8; zero new
    ... %              tunable parameters)
    ... %   quadMani - quadratic 1-D latent map for curved-curve modes
    ... %              (MRCoSEAManifoldFitQuad/SampleQuad; same 1-varThresh
    ... %              validity bar; m >= 2 modes keep incumbent behavior)
    def.segMani  = false;
    def.quadMani = false;

    cfg = def;
    if isstruct(userCfg)
        f = fieldnames(userCfg);
        for i = 1 : numel(f)
            cfg.(f{i}) = userCfg.(f{i});
        end
    end
end
