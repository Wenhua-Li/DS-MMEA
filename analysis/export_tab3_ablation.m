function export_tab3_ablation()
%EXPORT_TAB3_ABLATION Paper Tab.3: paired ablation of MSCoSEA (IGDX).
%   Two measured pairs, run-to-run seeded pairing on every problem:
%     - MSCoSEA vs MRCoSEA31 (previous generation, single-family reference);
%     - MSCoSEA vs MSCoSEA_nomani (manifold channel off).
%   For each problem the runs available in BOTH configs are paired (the 9
%   Polygon problems have 20 reinforced runs for MSCoSEA and MRCoSEA31;
%   the nomani arm only has 10 runs everywhere). Better/worse is decided by
%   the per-problem IGDX median; p comes from a Wilcoxon signed-rank test
%   on the paired runs of that problem.
%   Two further rows quote numbers recorded in the dev plan (Sec. 9, P2
%   pilot): the isolated-isle arm and the three-scale (K3) arm. Their raw
%   data was cleaned after the pilot; only the 24-problem per-run counts
%   (24 problems x 10 runs = 240 pairs) remain on record, which is a
%   different scope (per-run counts, not per-problem counts) and is marked
%   as such in the table.
%   Output: ../tab3_ablation.md
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); export_tab3_ablation"

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);       % repository root
    outPath = fullfile(here, 'out', 'tab3_ablation.md');
    if exist(fullfile(here,'out'),'dir') ~= 7, mkdir(fullfile(here,'out')); end
    dataRoot = fullfile(root, 'results','Data');

    probList = rq3_problems();
    probs  = {probList.name};
    groups = {probList.group};
    nP = numel(probs);

    pairs = {'MRCoSEA31', 44, 27; 'MSCoSEA_nomani', 46, 25};

    res = struct('name',{},'w',{},'t',{},'l',{},'pMed',{},'nSig',{}, ...
                 'nProb',{},'nPairs',{},'w10',{},'l10',{});
    for k = 1 : size(pairs, 1)
        w = 0; t = 0; l = 0; pv = nan(nP,1); nR = zeros(nP,1);
        for ip = 1 : nP
            a = localLoad(dataRoot, 'MSCoSEA', probs{ip}, 20);
            b = localLoad(dataRoot, pairs{k,1}, probs{ip}, 20);
            nCommon = min(numel(a), numel(b));
            nR(ip) = nCommon;
            a = a(1:nCommon); b = b(1:nCommon);
            if nCommon == 0, continue; end
            ma = median(a); mb = median(b);
            if ma < mb, w = w + 1;
            elseif ma > mb, l = l + 1;
            else, t = t + 1;
            end
            try
                pv(ip) = signrank(a, b);
            catch
                pv(ip) = 1;
            end
        end
        res(k) = struct('name',pairs{k,1},'w',w,'t',t,'l',l, ...
            'pMed',median(pv(isfinite(pv))),'nSig',sum(pv < 0.05), ...
            'nProb',sum(nR > 0),'nPairs',sum(nR),'w10',0,'l10',0);
        fprintf(['pair %d: %-15s better=%d tie=%d worse=%d  pMed=%.3g  ' ...
                 'sig=%d  maxRuns=%d\n'], k, pairs{k,1}, w, t, l, ...
                res(k).pMed, res(k).nSig, max(nR));
        % cross-check with the plain 10-run scope on all 71 problems
        w10 = 0; l10 = 0;
        for ip = 1 : nP
            a = localLoad(dataRoot, 'MSCoSEA', probs{ip}, 10);
            b = localLoad(dataRoot, pairs{k,1}, probs{ip}, 10);
            if isempty(a) || isempty(b), continue; end
            ma = median(a); mb = median(b);
            if ma < mb, w10 = w10 + 1;
            elseif ma > mb, l10 = l10 + 1; end
        end
        res(k).w10 = w10; res(k).l10 = l10;
        fprintf('         cross-check, 10 runs on all 71: better=%d worse=%d\n', ...
            w10, l10);
    end

    %% stdout sanity: expected anchors
    fprintf('\nanchor check: MRCoSEA31 expect 44 better / 27 worse; ');
    fprintf('nomani expect 46 better / 25 worse (v4.2 same-version, 2026-10-01)\n');

    %% markdown table
    L = {};
    L = [L; '# Tab.3  Paired ablation of MSCoSEA (IGDX, per-run seeded pairing)']; L = [L; ''];
    L = [L; '| Comparison | Scope | Better | Tie | Worse | median p | #(p<0.05) | Source |'];
    L = [L; '|---|---|---|---|---|---|---|---|'];
    for k = 1 : numel(res)
        L = [L; sprintf('| MSCoSEA vs `%s` | 71 problems | **%d** | %d | **%d** | %s | %d | measured (frozen runs) |', ...
             res(k).name, res(k).w, res(k).t, res(k).l, ...
             localFmtP(res(k).pMed), res(k).nSig)];
    end
    L = [L; '| MSCoSEA vs isle arm (no cross-family offspring) | 24-problem P2 pilot | **123** | — | **117** | — | — | dev plan §9 (data cleaned) |'];
    L = [L; '| MSCoSEA vs three-scale arm (K3) | 24-problem P2 pilot | **124** | — | **116** | — | — | dev plan §9 (data cleaned) |'];
    L = [L; '']; L = [L; '**Notes.**']; L = [L; ''];
    L = [L; ['1. Measured rows: for every problem the runs shared by both ' ...
             'configurations are compared pairwise (identical seeds). ' ...
             'Better/worse is decided on the per-problem IGDX medians; the ' ...
             'p-value is a per-problem Wilcoxon signed-rank test on the ' ...
             'paired runs (median over problems reported; the test returns ' ...
             'p = 1 when all paired differences are zero).']];
    L = [L; ['2. Run counts: 10 runs per problem everywhere, except the 9 ' ...
             'Polygon problems, which have 20 reinforced runs for ' ...
             'MSCoSEA and MRCoSEA31 and use all 20 in the pairing. The ' ...
             'nomani arm has 10 runs per problem, so that pair is 10-run ' ...
             'throughout.']];
    L = [L; ['3. Isle and three-scale rows: the arm data was cleaned after ' ...
             'the 24-problem P2 pilot; the counts quoted here are the ' ...
             'per-run paired counts recorded in the development plan (240 ' ...
             'run pairs = 24 problems x 10 runs), a different scope from ' ...
             'the per-problem counts above, kept for the contribution-chain ' ...
             'record only. The three-scale row direction follows the ' ...
             'MSCoSEA perspective of the recorded "K3 arm 116/124".']];
    L = [L; ['4. Ties can only appear when two medians coincide exactly; ' ...
             'no ties occurred in the measured rows.']];
    L = [L; ''];
    L = [L; '---']; L = [L; ''];

    %% write
    fid = fopen(outPath, 'w', 'n', 'UTF-8');
    assert(fid > 0, 'cannot open output file');
    fwrite(fid, strjoin(L, newline), 'char');
    fclose(fid);
    fprintf('tab3: wrote %s\n', outPath);
end

function v = localLoad(dataRoot, cfgName, probName, maxRuns)
%localLoad Collect the IGDX of the available runs r1..maxRuns in run order.
    v = [];
    for r = 1 : maxRuns
        f = fullfile(dataRoot, cfgName, sprintf('%s__r%d.mat', probName, r));
        if exist(f, 'file') ~= 2, continue; end
        s = load(f, 'IGDX');
        v(end+1) = s.IGDX;                                 %#ok<AGROW>
    end
    v = v(:);
end

function s = localFmtP(p)
%localFmtP Format a p-value for the table.
    if isnan(p)
        s = '—';
    elseif p < 0.001
        s = '<0.001';
    else
        s = sprintf('%.3g', p);
    end
end
