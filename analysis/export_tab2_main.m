function export_tab2_main()
%EXPORT_TAB2_MAIN Paper Tab.2: main comparison, 71 problems x 8 algorithms.
%   Methodology is bit-compatible with the frozen 71-problem report
%   (frozen internal 71-problem report, 2026-09-29):
%     - per problem, the IGDX median over the 10 paired runs;
%     - Friedman average ranks per problem (ties share the average rank;
%       a DNF cell, i.e. all runs failed -> NaN, takes the worst value);
%     - overall rank = mean rank over all 71 problems; group rank = mean
%       rank within each group G1-G6;
%     - head-to-head against MSCoSEA: per-problem median comparison for
%       better/tie/worse, plus the per-problem Wilcoxon signed-rank test on
%       the paired runs; the reported p is the median of the 71 p-values
%       (a rival DNF counts as a win with p undefined -> 1, as in the
%       internal report).
%   The 8 main configs: MSCoSEA + HREA, MMOEADC, CMMO, TriMOEATAR,
%   HHCMMEA, HDMMODE, FPITSEA (no CoMMEA, no in-house older versions).
%   Output: ../tab2_main.md ; the group-rank matrix is printed to stdout
%   for registration in paper/tables/numbers_map.md.
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); export_tab2_main"

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);       % repository root
    outPath = fullfile(here, 'out', 'tab2_main.md');
    if exist(fullfile(here,'out'),'dir') ~= 7, mkdir(fullfile(here,'out')); end
    dataRoot = fullfile(root, 'results','Data');

    probList = rq3_problems();
    probs  = {probList.name};
    groups = {probList.group};
    nP = numel(probs);
    runs = 1:10;

    %% 10 configs in the frozen order; main comparison = indices [1 4:10]
    cfg = {'MSCoSEA','MSCoSEA_nomani','MRCoSEA31', ...
           'HREA','MMOEADC','CMMO','TriMOEATAR','HHCMMEA','HDMMODE','FPITSEA'};
    mainIdx = [1 4 5 6 7 8 9 10];
    dispName = {'MSCoSEA','HREA','MMOEADC','CMMO','TriMOEA-TA&R', ...
                'HHC-MMEA','HDMMODE','FPITSEA'};
    nM = numel(mainIdx);

    %% load paired runs (IGDX only)
    RAW = cell(nM, nP);
    MED = nan(nM, nP);
    for k = 1 : nM
        for ip = 1 : nP
            vx = [];
            for r = runs
                f = fullfile(dataRoot, cfg{mainIdx(k)}, ...
                             sprintf('%s__r%d.mat', probs{ip}, r));
                if exist(f, 'file') ~= 2, continue; end
                s = load(f, 'IGDX');
                vx(end+1) = s.IGDX;                        %#ok<AGROW>
            end
            if numel(vx) == 0, continue; end               % DNF cell
            RAW{k,ip} = vx(:);
            MED(k,ip) = median(vx);
        end
    end
    nMiss = sum(~isfinite(MED(:)));
    if nMiss > 0
        fprintf('note: %d DNF cells handled as worst rank\n', nMiss);
    end

    %% Friedman average ranks (identical routine to the internal report)
    R = localRanks(MED);
    ovm = mean(R, 2);
    gTags = {'G1','G2','G3','G4','G5','G6'};
    gRm = zeros(nM, numel(gTags));
    for gi = 1 : numel(gTags)
        gRm(:,gi) = mean(R(:, strcmp(groups, gTags{gi})), 2);
    end

    %% head-to-head vs MSCoSEA
    W = zeros(nM,1); T = zeros(nM,1); Lo = zeros(nM,1);
    pMed = nan(nM,1); nSig = nan(nM,1);
    for k = 1 : nM
        if k == 1, continue; end
        w = 0; t = 0; l = 0; pv = nan(nP,1);
        for ip = 1 : nP
            a = MED(1,ip); b = MED(k,ip);
            if ~isfinite(b)                    % rival DNF on this problem
                w = w + 1;
            elseif a < b
                w = w + 1;
            elseif a > b
                l = l + 1;
            else
                t = t + 1;
            end
            try
                pv(ip) = signrank(RAW{1,ip}, RAW{k,ip});
            catch
                pv(ip) = 1;
            end
        end
        W(k) = w; T(k) = t; Lo(k) = l;
        pMed(k) = median(pv);
        nSig(k) = sum(pv < 0.05);
    end

    %% stdout: group-rank matrix (registration copy)
    fprintf('\n== group-rank matrix (IGDX mean rank, lower is better) ==\n');
    fprintf('%-14s', 'algorithm');
    fprintf('%7s', gTags{:}); fprintf('%9s\n', 'Overall');
    [ ~, ord ] = sort(ovm);
    for k = ord.'                        % row vector: iterate scalars
        fprintf('%-14s', dispName{k});
        fprintf('%7.2f', gRm(k,:));
        fprintf('%9.2f\n', ovm(k));
    end
    fprintf('\nanchor check: MSCoSEA overall = %.2f (expect 2.72), ', ovm(1));
    iMM = find(strcmp(cfg(mainIdx), 'MMOEADC'));
    fprintf('MMOEADC overall = %.2f (expect 3.54)\n', ovm(iMM));

    %% markdown table
    L = {};
    L = [L; '# Tab.2  Main comparison on 71 problems (IGDX, median over 10 paired runs)']; L = [L; ''];
    L = [L; '| Algorithm | MMF | IDMP | Local-front | Polygon | Mode-count | HDMMF | Overall | W/T/L vs MSCoSEA | median p | #(p<0.05) |'];
    L = [L; '|---|---|---|---|---|---|---|---|---|---|---|'];
    for k = ord.'
        if ovm(k) <= min(ovm) + 1e-9, nm = ['**' dispName{k} '**']; else, nm = dispName{k}; end
        if k == 1
            h2h = '—'; ps = '—'; ns = '—';
        else
            h2h = sprintf('%d/%d/%d', W(k), T(k), Lo(k));
            ps = localFmtP(pMed(k));
            ns = sprintf('%d', nSig(k));
        end
        L = [L; sprintf('| %s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | **%.2f** | %s | %s | %s |', ...
             nm, gRm(k,:), ovm(k), h2h, ps, ns)];
    end
    L = [L; '']; L = [L; '**Notes.**']; L = [L; ''];
    L = [L; ['1. Metric: IGDX (decision-space IGD), median of the 10 ' ...
             'independent runs per problem; all configurations use the ' ...
             'identical seed table (20260907 + 1000·r + problem index).']];
    L = [L; ['2. Ranks: Friedman average ranks per problem (ties share the ' ...
             'average rank). DNF = all runs crashed (MMOEADC on NMMF5; ' ...
             'objective values overflow deterministically in its grid-based ' ...
             'environmental selection); the cell is placed last (worst rank) ' ...
             'and counts as a win for MSCoSEA in the head-to-head.']];
    L = [L; ['3. W/T/L = number of problems where MSCoSEA''s median IGDX ' ...
             'is lower / equal / higher than the compared algorithm''s ' ...
             '(i.e., wins / ties / losses from the MSCoSEA perspective). ' ...
             'The p column is the median over the 71 per-problem Wilcoxon ' ...
             'signed-rank tests on the paired runs; #(p<0.05) counts the ' ...
             'problems with a significant test.']];
    L = [L; ['4. Group sizes: MMF = 13, IDMP = 12, Local-front = 16, ' ...
             'Polygon = 9, Mode-count = 6, HDMMF = 15. Overall = mean rank ' ...
             'over all 71 problems.']];
    L = [L; ''];
    L = [L; '---']; L = [L; ''];

    %% write
    fid = fopen(outPath, 'w', 'n', 'UTF-8');
    assert(fid > 0, 'cannot open output file');
    fwrite(fid, strjoin(L, newline), 'char');
    fclose(fid);
    fprintf('tab2: wrote %s\n', outPath);
end

function R = localRanks(M)
%localRanks Friedman-style average ranks per column (small = good).
%   Non-finite entries (NaN = DNF) are treated as the worst value.
    R = zeros(size(M));
    for ip = 1 : size(M, 2)
        tbl = M(:, ip);
        tbl(~isfinite(tbl)) = inf;
        [~, rk] = sort(tbl);
        [~, rk] = sort(rk);
        for vv = unique(tbl(:)).'
            idx = tbl == vv;
            if sum(idx) > 1
                rk(idx) = mean(rk(idx));
            end
        end
        R(:, ip) = rk;
    end
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
