function retain = CoSEAEpsRetain(PopObj, FrontNo, epsT)
%CoSEAEpsRetain Local-PS retention via multiplicative epsilon-dominance (M8).
%   retain = CoSEAEpsRetain(PopObj, FrontNo, epsT) marks individuals outside
%   the first front that are NOT dominated by the scaled first front
%   (1+epsT)*F1 (multiplicative epsilon-dominance, CoMMEA/CMMO style). Such
%   individuals are approximately objective-equivalent to the front but sit
%   elsewhere in the decision space, i.e. candidate newborn/rare modes.
%   Assumes minimization with non-negative objectives (holds for MMF/IDMP).

    n      = size(PopObj, 1);
    retain = false(1, n);
    F1     = PopObj(FrontNo == 1, :);
    if isempty(F1) || epsT <= 0
        return;
    end
    scaled  = F1 * (1 + epsT);
    candIdx = find(FrontNo > 1);
    for i = candIdx(:).'
        x   = PopObj(i, :);
        dom = all(scaled <= x, 2) & any(scaled < x, 2);
        if ~any(dom)
            retain(i) = true;
        end
    end
end
