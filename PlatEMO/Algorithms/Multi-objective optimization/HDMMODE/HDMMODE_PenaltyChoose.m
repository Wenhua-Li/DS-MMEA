function Penalty_index = HDMMODE_PenaltyChoose(Pop,Archive,rr)
%HDMMODE_PenaltyChoose Indices of population members covered by the archive.
%   A member is "covered" when an archive point is within rr(i) of it in
%   every dimension; such members are penalized in HDMMODE_EnvSelect2.
%
%   PORT FIX (2026-09-27): the original code used the linear indices of
%   find(D == dim) directly as ROW indices of the population. Whenever a
%   match sat in a column j > 1 the linear index exceeded the number of
%   rows and the call crashed (or silently wrapped). Fixed to the intended
%   row semantics: unique row indices of the matching entries.

    Pop_dec    = Pop.decs;
    Archive_dec = Archive.decs;
    [n1,dim] = size(Pop_dec);
    n2 = size(Archive_dec,1);
    D = zeros(n1,n2);
    for i = 1 : dim
        D = D + (abs(Pop_dec(:,i) - repmat(Archive_dec(:,i)',n1,1)) <= rr(i));
    end
    [row,~] = find(D == dim);
    Penalty_index = unique(row);
end
