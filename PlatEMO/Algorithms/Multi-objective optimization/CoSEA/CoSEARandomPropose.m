function cand = CoSEARandomPropose(Problem,nCand)
%CoSEARandomPropose Uniform random proposals in the decision box (A3 control).
%   cand = CoSEARandomPropose(Problem,nCand) returns nCand x D points drawn
%   uniformly from the decision box. Direction-free control for the
%   reflection proposal: identical injection schedule, identical number of
%   evaluated candidates (hence identical evaluation budget), but no
%   structural prior about where undiscovered equivalent Pareto sets may lie.

    lb = Problem.lower;
    ub = Problem.upper;
    cand = unifrnd(repmat(lb,nCand,1),repmat(ub,nCand,1));
end
