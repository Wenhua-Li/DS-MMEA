function Del = HDMMODE_Truncation2(PopObj,K)
%HDMMODE_Truncation2 Select K solutions for deletion by distance truncation.
%   Iteratively removes the solution whose sorted distance vector is
%   lexicographically smallest (SPEA2-style).
%   Ported from the original HDMMODE code (renamed); logic unchanged.

    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain = find(~Del);
        Temp   = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
