function modes = CoSEATruePSModes(Problem, N)
%CoSEATruePSModes True Pareto sets per mode (decision space) of a problem.
%   modes = CoSEATruePSModes(Problem, N) returns a cell array of per-mode
%   true-PS samples. Uses Problem.GetPSModes(N) when available (IDMP_p);
%   otherwise clusters Problem.GetOptimum(N) into modes by connectivity in
%   the box-normalized decision space (radius 0.05*sqrt(D)).

    if nargin < 2
        N = 2000;
    end
    if ismethod(Problem, 'GetPSModes')
        modes = Problem.GetPSModes(N);
        return;
    end
    if isprop(Problem, 'POS')
        PS = Problem.POS;               % cached true Pareto set (MMF family)
    else
        PS = Problem.optimum;           % assume decision-space optima
    end
    if size(PS, 1) > N                  % subsample for speed
        PS = PS(round(linspace(1, size(PS,1), N)), :);
    end
    lb = Problem.lower;
    ub = Problem.upper;
    PSn = (PS - lb) ./ max(ub - lb, eps);
    thr = 0.05 * sqrt(Problem.D);
    Dm  = pdist2(PSn, PSn);
    adj = Dm < thr;
    n   = size(PSn, 1);
    adj(1:n+1:end) = false;
    comp  = conncomp(graph(sparse(adj)));
    nComp = max(comp);
    modes = cell(1, nComp);
    for c = 1 : nComp
        modes{c} = PS(comp == c, :);
    end
end
