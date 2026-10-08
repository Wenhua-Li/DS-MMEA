function T = CoSEATruePS(Problem,probName,N,cacheDir)
%CoSEATruePS Unified ground truth (true PS, true PF, per-mode PS) of a problem.
%   T = CoSEATruePS(Problem,probName) returns a struct with
%     T.PS     - reference true Pareto set  (n x D, raw decision space)
%     T.PF     - reference true Pareto front (m x M, raw objective space)
%     T.modes  - cell array, one cell per equivalent Pareto set (mode)
%     T.nModes - number of modes
%     T.src    - where the truth came from (audit trail)
%     T.D/T.M/T.lower/T.upper/T.name
%   T = CoSEATruePS(Problem,probName,N) sets the sampling size for analytic
%   mode generators (default 2000). T = CoSEATruePS(...,cacheDir) overrides
%   the on-disk cache directory.
%
%   Source priority (fixed so that all algorithms are scored against the very
%   same reference sets, and so that numbers stay comparable with the
%   published CoMMEA/HREA protocol):
%     1. Problem.GetPSModes(N)  - exact analytic modes (IDMP_p, Polygon_p)
%     2. reference data files   - PlatEMO/ReferenceData/, the published PS/PF
%        samples: Reference_PSPF_data/<name>_Reference_PSPF_data.mat (fields
%        PS,PF), TruePSPFdata/<name>.mat (IDMP: PSS,PF; IDMP_e: PSs,PFs; MMF:
%        PS,PF via <name>truePSPF.mat), PFPS/<name>.mat (PSS,PF)
%     3. Problem.POS            - analytic PS cached by the problem class
%   When the source does not provide modes (cases 2 and 3) the PS sample is
%   split into modes by connectivity clustering (CoSEAClusterModes).
%
%   Results are cached as <cacheDir>/<probName>__N<N>.mat. Build the cache
%   serially first (exp/rq3/build_truth_cache.m); parallel workers then only
%   read it. Writes go to a temp file and are moved into place, so a partial
%   file can never be read by another worker.

    if nargin < 3 || isempty(N)
        N = 2000;
    end
    here = fileparts(mfilename('fullpath'));
    if nargin < 4 || isempty(cacheDir)
        cacheDir = fullfile(fileparts(here),'exp','rq3','TruthCache');
    end
    if exist(cacheDir,'dir') ~= 7
        mkdir(cacheDir);
    end
    cacheFile = fullfile(cacheDir,sprintf('%s__N%d.mat',probName,N));
    if exist(cacheFile,'file') == 2
        try
            S = load(cacheFile);
            if isfield(S,'T') && isfield(S.T,'PS') && isfield(S.T,'PF') && ...
                    isfield(S.T,'modes') && ~isempty(S.T.PS)
                T = S.T;
                return;
            end
        catch
        end
    end

    T = localBuild(Problem,probName,N,here);
    try
        tmp = sprintf('%s.%d.tmp',cacheFile,feature('getpid'));
        save(tmp,'T','-v7');
        movefile(tmp,cacheFile,'f');
    catch
    end
end

function T = localBuild(Problem,probName,N,here)
    refRoot = fullfile(fileparts(here),'PlatEMO','ReferenceData');
    T.name   = probName;
    T.D      = Problem.D;
    T.M      = Problem.M;
    T.lower  = Problem.lower;
    T.upper  = Problem.upper;

    modes = {};
    if ismethod(Problem,'GetPSModes')
        modes = Problem.GetPSModes(N);
        modes = modes(~cellfun(@isempty,modes));
        T.PS  = vertcat(modes{:});
        T.PF  = Problem.PF;
        T.src = 'GetPSModes';
    else
        [PS,PF,src] = localLoadReference(refRoot,probName,T.D,T.M);
        if isempty(PS) && isprop(Problem,'POS') && ~isempty(Problem.POS)
            PS  = Problem.POS;
            PF  = Problem.PF;
            src = 'POS';
        end
        if isempty(PS)
            error('CoSEATruePS:noTruth', ...
                'No ground truth available for problem "%s".',probName);
        end
        T.PS  = PS;
        T.PF  = PF;
        T.src = src;
    end

    if isempty(modes)
        modes = CoSEAClusterModes(T.PS,T.lower,T.upper,0.05*sqrt(T.D));
    end
    T.modes  = modes;
    T.nModes = numel(modes);

    % Mode resolvability audit: the coverage indicators (nicheRecovery, the
    % normalized rPSP) declare a mode recovered when the population comes
    % within delta = 0.05 of it in the BOX-normalized decision space. When the
    % true PS occupies only a small part of the decision box (the polygon
    % group: box [-100,100]^D, PS extent ~12), distinct modes can be closer
    % than delta and the recovery indicator degenerates to 1.0 for every
    % algorithm. T.dSep is the smallest normalized distance between two
    % distinct modes; T.resolvable tells the evaluation layer whether the
    % recovery indicators are meaningful for this problem (raw-space IGD/IGDX/
    % CR/PSP are unaffected and always reported).
    T.delta = 0.05;
    T.dSep  = localModeSep(modes,T.lower,T.upper);
    T.resolvable = (T.nModes < 2) || (T.dSep > T.delta);
end

function d = localModeSep(modes,lower,upper)
%localModeSep Smallest box-normalized distance between two distinct modes.
    K = numel(modes);
    if K < 2
        d = Inf;
        return;
    end
    sp = max(upper - lower,eps);
    Mn = cell(1,K);
    for j = 1 : K
        X = (modes{j} - lower) ./ sp;
        if size(X,1) > 300
            idx = unique(round(linspace(1,size(X,1),300)));
            X   = X(idx,:);
        end
        Mn{j} = X;
    end
    d = Inf;
    for a = 1 : K-1
        for b = a+1 : K
            d = min(d,min(pdist2(Mn{a},Mn{b}),[],'all'));
        end
    end
end

function [PS,PF,src] = localLoadReference(refRoot,probName,D,M)
    PS = []; PF = []; src = '';
    cands = { ...
        fullfile(refRoot,'Reference_PSPF_data',[probName '_Reference_PSPF_data.mat']); ...
        fullfile(refRoot,'TruePSPFdata',[probName '.mat']); ...
        fullfile(refRoot,'TruePSPFdata',[probName 'truePSPF.mat']); ...
        fullfile(refRoot,'PFPS',[probName '.mat'])};
    psFields = {'PSS','PSs','PS'};      % dense sample first (IDMP: PS is 2x2)
    pfFields = {'PF','PFs'};
    for i = 1 : size(cands,1)
        f = cands{i,1};
        if exist(f,'file') ~= 2
            continue;
        end
        S = load(f);
        ps = []; pf = [];
        for k = 1 : numel(psFields)
            if isfield(S,psFields{k}) && ~isempty(S.(psFields{k}))
                ps = S.(psFields{k});
                break;
            end
        end
        for k = 1 : numel(pfFields)
            if isfield(S,pfFields{k}) && ~isempty(S.(pfFields{k}))
                pf = S.(pfFields{k});
                break;
            end
        end
        if isempty(ps) || isempty(pf)
            continue;
        end
        if size(ps,2) ~= D || size(pf,2) ~= M
            continue;                    % dimension mismatch: not this problem
        end
        PS  = ps;
        PF  = pf;
        rel = strrep(f,[refRoot filesep],'');
        src = sprintf('file:%s',strrep(rel,filesep,'/'));
        return;
    end
end
