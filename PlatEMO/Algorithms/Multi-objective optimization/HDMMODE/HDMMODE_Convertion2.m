function outcome = HDMMODE_Convertion2(Average,G,last_gen,change_threshold)
%HDMMODE_Convertion2 Stagnation test on the normalized-objective average.
%   Returns true when the recorded average has changed by at most
%   change_threshold over the last last_gen generations.
%   Ported from the original HDMMODE code (renamed); logic unchanged.

    if max(abs(Average(G,:) - Average(G-last_gen,:))) <= change_threshold
        outcome = true;
    else
        outcome = false;
    end
end
