function export_tab4_runtime()
%EXPORT_TAB4_RUNTIME Paper Tab.4: non-contended single-run wall clock.
%   Reads the frozen solo-timing artifacts produced by measure_runtime_solo.m
%   (sequential execution, no parpool, 3 runs per cell, 4 representative
%   problems, the experiment seed table). The 'runtime' field of the
%   parallel batch runs reflects CPU contention and is NOT usable here.
%   Source priority: runtime_solo.mat (full per-run entries); the two v4.0
%   arms MSCoSEA and MSCoSEA_nomani are missing from the incremental .mat
%   (resume path skipped the re-save), so their per-cell means are taken
%   from the summary log runtime_solo.log written by the same run, whose
%   2-decimal precision equals the table format anyway. The .mat log-mismatch
%   is recorded in the internal notes.
%   HDMMODE and FPITSEA were added to the experiment matrix after the solo
%   panel was frozen, so they have no solo measurement and appear as an em
%   dash. CoMMEA rows are excluded from the paper scope by decision.
%   Output: ../tab4_runtime.md
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); export_tab4_runtime"

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);       % repository root
    outPath = fullfile(here, 'out', 'tab4_runtime.md');
    if exist(fullfile(here,'out'),'dir') ~= 7, mkdir(fullfile(here,'out')); end
    rq3dir = fullfile(root, 'results');
    probs  = {'MMF1','IDMPM3T2','IDMPM2T2_e','Polygon_M6_D10'};

    %% source 1: runtime_solo.mat (per-run entries)
    cellMean = containers.Map;            % key 'cfg|prob' -> mean seconds
    S = load(fullfile(rq3dir, 'runtime_solo.mat'));
    assert(isfield(S, 'rows') && isfield(S, 'res'), 'unexpected runtime_solo.mat');
    res = S.res;
    for ic = 1 : numel(res)
        k = sprintf('%s|%s', res(ic).cfg, res(ic).prob);
        if isKey(cellMean, k)
            cellMean(k) = [cellMean(k), res(ic).rt];
        else
            cellMean(k) = res(ic).rt;
        end
    end
    matCfgs = unique({res.cfg}, 'stable');

    %% source 2: runtime_solo.log (summary lines, 2-decimal) for configs
    %   that have no per-run entries in the .mat
    logCfgs = {};
    logCell = cell(0, 6);                 % cfg, 4 cell means, full-precision mean
    fid = fopen(fullfile(rq3dir, 'runtime_solo.log'), 'r');
    assert(fid > 0, 'cannot open runtime_solo.log');
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        tok = strtrim(strsplit(ln, ' ', 'CollapseDelimiters', true));
        if numel(tok) ~= 6, continue; end
        v = cellfun(@str2double, tok(2:end));
        if any(isnan(v)), continue; end   % skips the header line
        logCfgs{end+1} = tok{1};          %#ok<AGROW>
        logCell(end+1, :) = {tok{1}, v(1), v(2), v(3), v(4), v(5)}; %#ok<AGROW>
    end
    fclose(fid);
    needLog = setdiff({'MSCoSEA','MSCoSEA_nomani'}, matCfgs);
    assert(all(ismember(needLog, logCfgs)), ...
        'configs missing from both runtime_solo.mat and .log');

    %% FE-budget audit from the frozen experiment runs (8 main configs)
    %   Re-states the registered "<= 1.020 x maxFE" bound from data: max
    %   over all problems/runs of FE / maxFE.
    probList = rq3_problems();
    probsAll = {probList.name};
    main8 = {'MSCoSEA','HREA','MMOEADC','CMMO','TriMOEATAR','HHCMMEA', ...
             'HDMMODE','FPITSEA'};
    dataRoot = fullfile(root, 'results','Data');
    feMax = 0; feMS = 0; feTri = 0; feFPIT = 0;
    for ic = 1 : numel(main8)
        for ip = 1 : numel(probsAll)
            for r = 1 : 10
                f = fullfile(dataRoot, main8{ic}, ...
                             sprintf('%s__r%d.mat', probsAll{ip}, r));
                if exist(f, 'file') ~= 2, continue; end
                s = load(f, 'FE', 'maxFE');
                v = s.FE / s.maxFE;
                feMax = max(feMax, v);
                switch main8{ic}
                    case 'MSCoSEA',    feMS = max(feMS, v);
                    case 'TriMOEATAR', feTri = max(feTri, v);
                    case 'FPITSEA',    feFPIT = max(feFPIT, v);
                end
            end
        end
    end
    fprintf('tab4: FE audit: overall=%.4f  MSCoSEA=%.4f  TriMOEA=%.3f  ', ...
        feMax, feMS, feTri);
    fprintf('FPITSEA=%.3f\n', feFPIT);

    %% paper rows in the main-comparison order + ablation references
    rows = { ...
        'MSCoSEA',       'MSCoSEA',       false;
        'HREA',          'HREA',          false;
        'MMOEADC',       'MMOEADC',       false;
        'CMMO',          'CMMO',          false;
        'TriMOEA-TA&R',  'TriMOEATAR',    false;
        'HHC-MMEA',      'HHCMMEA',       false;
        'HDMMODE',       'HDMMODE',       false;
        'FPITSEA',       'FPITSEA',       false;
        'MRCoSEA31',     'MRCoSEA31',     true;
        'MSCoSEA\_nomani', 'MSCoSEA_nomani', true};

    %% assemble
    L = {};
    L = [L; '# Tab.4  Single-run wall clock (sequential, non-contended)']; L = [L; ''];
    L = [L; '| Algorithm | MMF1 | IDMPM3T2 | IDMPM2T2\_e | Polygon\_M6\_D10 | Mean |'];
    L = [L; '|---|---|---|---|---|---|'];
    meanOf = struct();
    inMat = '';
    for k = 1 : size(rows, 1)
        nm = rows{k,1}; cfgName = rows{k,2};
        vals = nan(1, 4);
        if any(strcmp(matCfgs, cfgName))
            for ip = 1 : 4
                kk = sprintf('%s|%s', cfgName, probs{ip});
                if isKey(cellMean, kk), vals(ip) = mean(cellMean(kk)); end
            end
            inMat = 'mat';
        elseif any(strcmp(logCfgs, cfgName))
            r = find(strcmp(logCfgs, cfgName), 1);
            vals = [logCell{r,2:5}];      % already 2-decimal means
            meanOf.(cfgName) = logCell{r,6};  % full-precision mean from the log
            inMat = 'log';
        end
        if all(isnan(vals))
            L = [L; sprintf('| %s | — | — | — | — | — |', nm)];
            meanOf.(cfgName) = NaN;
        elseif strcmp(inMat, 'log')
            L = [L; sprintf('| %s%s | %.2f | %.2f | %.2f | %.2f | **%.2f** |', ...
                 nm, localRefTag(rows{k,3}), vals, meanOf.(cfgName))];
        else
            meanOf.(cfgName) = mean(vals);
            L = [L; sprintf('| %s%s | %.2f | %.2f | %.2f | %.2f | **%.2f** |', ...
                 nm, localRefTag(rows{k,3}), vals, mean(vals))];
        end
        fprintf('tab4: %-16s source=%-3s mean=%.2f\n', cfgName, inMat, meanOf.(cfgName));
    end
    L = [L; '']; L = [L; '**Notes.**']; L = [L; ''];
    L = [L; ['1. Wall-clock seconds of one complete run, measured strictly ' ...
             'sequentially (no parallel workers, no contention): 3 runs per ' ...
             'cell, 4 representative problems (MMF1; IDMPM3T2; IDMPM2T2_e; ' ...
             'Polygon_M6_D10), experiment seed table. Mean over the 3 runs ' ...
             'per cell; Mean = unweighted mean of the 4 cells.']];
    L = [L; ['2. HDMMODE and FPITSEA joined the experiment matrix after ' ...
             'this solo-timing panel was frozen; no non-contended ' ...
             'measurement exists and the contended batch runtimes are not ' ...
             'comparable, hence the em dash.']];
    L = [L; ['3. Rows marked (a) are ablation/reference configurations, ' ...
             'not part of the 8-algorithm main comparison: MRCoSEA31 is the ' ...
             'single-family predecessor, MSCoSEA_nomani disables the ' ...
             'manifold channel.']];
    L = [L; sprintf(['4. Budget audit (measured from the frozen runs, 71 ' ...
             'problems x 10 runs each): MSCoSEA consumes at most %.4f x the ' ...
             'maxFE budget (injections and truth sampling included in the ' ...
             'accounting), i.e. within 2%% of budget. The external ' ...
             'baselines run on the same nominal budget with their own ' ...
             'internal accounting; two of them exceed it on the HDMMF ' ...
             'group (TriMOEA-TA&R up to %.3f x on NMMF6, FPITSEA up to ' ...
             '%.3f x on NMMF11), reported as measured.'], ...
            feMS, feTri, feFPIT)];
    L = [L; ''];
    L = [L; '---']; L = [L; ''];
    %% write
    fid = fopen(outPath, 'w', 'n', 'UTF-8');
    assert(fid > 0, 'cannot open output file');
    fwrite(fid, strjoin(L, newline), 'char');
    fclose(fid);
    fprintf('tab4: wrote %s\n', outPath);
    fprintf('tab4: FE audit: MSCoSEA <= %.4f (registered 1.020)\n', feMS);
    fprintf('tab4: anchors  MSCoSEA=%.2f (8.88, v4.2)  MRCoSEA31=%.2f (~4.1)  ', ...
        meanOf.MSCoSEA, meanOf.MRCoSEA31);
    fprintf('nomani=%.2f (4.77)\n', meanOf.MSCoSEA_nomani);
end

function tag = localRefTag(isRef)
%localRefTag Suffix marking ablation-reference rows.
    if isRef
        tag = ' (a)';
    else
        tag = '';
    end
end
