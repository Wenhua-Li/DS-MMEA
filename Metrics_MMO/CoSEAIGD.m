function val = CoSEAIGD(Obtained,Reference)
%CoSEAIGD Inverted generational distance in the raw (unnormalized) space.
%   val = CoSEAIGD(Obtained,Reference) returns the mean over reference points
%   of the Euclidean distance to the nearest obtained point:
%     val = (1/|R|) * sum_{y in R} min_{x in Obtained} ||x - y||
%   Used for IGD (objective space) and IGDX (decision space). Identical to
%   IGD_calculation.m of the CoMMEA reference implementation (JAS 2023),
%   verified against its stored results in tests/test_metrics_std.m.
%
%   Obtained  - n x d matrix (obtained set)
%   Reference - m x d matrix (reference set, e.g. true PS or true PF)

    if isempty(Obtained)
        val = Inf;
        return;
    end
    Dm  = pdist2(Reference,Obtained);   % m x n
    val = mean(min(Dm,[],2));
end
