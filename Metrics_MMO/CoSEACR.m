function CR = CoSEACR(Obtained,Reference)
%CoSEACR Cover rate of an obtained set with respect to a reference set.
%   CR = CoSEACR(Obtained,Reference) is the geometric mean, over decision
%   space dimensions, of the interval overlap ratio between the bounding box
%   of Obtained and the bounding box of Reference:
%     kesi_i = ((min(omax_i,rmax_i) - max(omin_i,rmin_i)) / (rmax_i-rmin_i))^2
%     CR     = nthroot(prod(kesi), 2*n_var)          % = geomean of ratios
%   Dimensions in which the reference set is degenerate (rmax == rmin) are
%   ignored (kesi = 1); disjoint intervals give kesi = 0 hence CR = 0.
%   Faithful port of CR_calculation.m of the CoMMEA reference implementation
%   (JAS 2023); verified in tests/test_metrics_std.m. PSP = CR/IGDX (Eq. 8 of
%   that paper) is assembled in CoSEAMetricsStd.
%
%   Obtained  - n x d obtained decision vectors
%   Reference - m x d reference (true) decision vectors

    n_var = size(Reference,2);
    if size(Obtained,1) == 1
        omin = Obtained;
        omax = Obtained;
    else
        omin = min(Obtained,[],1);
        omax = max(Obtained,[],1);
    end
    rmin = min(Reference,[],1);
    rmax = max(Reference,[],1);

    kesi = ones(1,n_var);
    for i = 1 : n_var
        if rmax(i) == rmin(i)
            kesi(i) = 1;                       % degenerate reference dimension
        elseif omin(i) >= rmax(i) || rmin(i) >= omax(i)
            kesi(i) = 0;                       % disjoint intervals
        else
            kesi(i) = ((min(omax(i),rmax(i)) - max(omin(i),rmin(i))) / ...
                       (rmax(i) - rmin(i)))^2;
        end
    end
    kesi = min(max(kesi,0),1);                 % guard round-off
    CR   = nthroot(prod(kesi),2*n_var);
end
