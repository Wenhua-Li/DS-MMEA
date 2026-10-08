function S = CoSEAMetricsStd(PopDec,PopObj,T)
%CoSEAMetricsStd Literature-standard MMO indicators in the raw space.
%   S = CoSEAMetricsStd(PopDec,PopObj,T) returns
%     S.IGD   - inverted generational distance in the objective space
%               (obtained objectives vs reference PF T.PF)
%     S.IGDX  - inverted generational distance in the decision space
%               (obtained decisions vs reference PS T.PS)
%     S.CR    - cover rate (CoSEACR) in the decision space
%     S.PSP   - CR / IGDX   (large is better, Eq. 8 of CoMMEA JAS 2023)
%     S.rPSP  - IGDX / CR   (small is better; Inf when CR = 0)
%   T is the ground-truth struct produced by CoSEATruePS (fields PS, PF).
%   All quantities are computed in the raw (unnormalized) space so that the
%   numbers are directly comparable with published CoMMEA/HREA tables. The
%   box-normalized coverage indicators (rPSP.m, nicheRecovery.m) are computed
%   separately by the pipeline.

    S.IGDX = CoSEAIGD(PopDec,T.PS);
    S.IGD  = CoSEAIGD(PopObj,T.PF);
    S.CR   = CoSEACR(PopDec,T.PS);
    if S.CR > 0 && S.IGDX > 0 && isfinite(S.IGDX)
        S.PSP  = S.CR / S.IGDX;
        S.rPSP = S.IGDX / S.CR;
    else
        S.PSP  = 0;
        S.rPSP = Inf;
    end
end
