function [recRate, recFlag, fullRec] = nicheRecovery(PopDec, TruePSModes, lower, upper, delta)
%nicheRecovery Niche (equivalent Pareto set) recovery rate.
%   [recRate,recFlag,fullRec] = nicheRecovery(PopDec, TruePSModes, lower, upper, delta)
%   TruePSModes is a cell array whose j-th cell contains sampled points of
%   the j-th true Pareto set (decision space). Mode j is recovered when its
%   nearest true-PS point is within delta (box-normalized Euclidean) of the
%   obtained population. recRate = mean(recFlag); fullRec = all(recFlag).

    if nargin < 5 || isempty(delta)
        delta = 0.05;
    end
    K = numel(TruePSModes);
    sp = max(upper - lower, eps);
    Pn = (PopDec - lower) ./ sp;
    recFlag = false(1, K);
    for j = 1 : K
        Tn = (TruePSModes{j} - lower) ./ sp;
        dMin = min(pdist2(Tn, Pn), [], 2);
        recFlag(j) = min(dMin) <= delta;
    end
    recRate = mean(recFlag);
    fullRec = all(recFlag);
end
