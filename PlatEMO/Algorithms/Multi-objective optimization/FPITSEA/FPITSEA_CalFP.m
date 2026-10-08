function IFP = FPITSEA_CalFP(Union)
% Calculate the fuzzy preference indicator for each solution in Union.
% Renamed from CalFP to avoid potential conflicts across algorithms.

    IFP = zeros(1,length(Union));

    %% Convergence indicator (count of dominators)
    Obj = Union.objs;
    n   = length(Union);
    Dominate = false(n);
    for i = 1 : n-1
        for j = i+1 : n
            k = any(Obj(i,:)<Obj(j,:)) - any(Obj(i,:)>Obj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    Conv = sum(Dominate);
    Con  = (Conv - min(Conv)) ./ (max(Conv) - min(Conv));

    %% Diversity indicator (decision-space crowding)
    CD  = FPITSEA_Crowding(Union.decs);
    Div = (CD - min(CD)) ./ (max(CD) - min(CD));

    %% Fuzzy preference indicator
    Con = Con';
    Div = Div';
    for i = 1 : n
        A = 1:n;
        A = setdiff(A,i);
        d_div12 = (Div(i) - Div(A))';
        d_con12 = (Con(A) - Con(i))';
        udiv12  = 1 ./ (1 + exp(-d_div12./1));
        ucon12  = 1 ./ (1 + exp(-d_con12./1));
        tmp = sum(udiv12 .* ucon12);
        IFP(i) = tmp;
    end
end