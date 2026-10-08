function [rPSPval, PSPval, CRval, IGDXval] = rPSP(PopDec, TruePS, lower, upper, delta)
%rPSP Reciprocal Pareto Set Proximity (small is better) in the decision space.
%   [rPSP,PSP,CR,IGDX] = rPSP(PopDec, TruePS, lower, upper, delta) computes
%     IGDX = mean over true-PS points of the distance to the nearest PopDec
%     CR   = fraction of PopDec points within delta of the true PS
%     PSP  = CR / IGDX   (large is better)
%     rPSP = IGDX / CR   (small is better; Inf when CR = 0)
%   All distances are Euclidean in the box-normalized decision space.
%   Note: this is the self-consistent definition of the rebuilt pipeline;
%   the historical CODES-thesis implementation was lost, so absolute values
%   are only comparable within the new pipeline, not to old logs.

    if nargin < 5 || isempty(delta)
        delta = 0.05;
    end
    sp   = max(upper - lower, eps);
    Pn   = (PopDec - lower) ./ sp;
    Tn   = (TruePS - lower) ./ sp;
    Dm   = pdist2(Tn, Pn);                    % |PS| x |Pop|
    dMin = min(Dm, [], 2);                    % PS side: nearest sample
    IGDXval = mean(dMin);
    CRval   = mean(min(Dm, [], 1) <= delta);  % Pop side: near-PS fraction
    if CRval > 0
        PSPval  = CRval / IGDXval;
        rPSPval = IGDXval / CRval;
    else
        PSPval  = 0;
        rPSPval = inf;
    end
end
