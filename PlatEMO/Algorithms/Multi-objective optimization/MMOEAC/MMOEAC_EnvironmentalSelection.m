function Population = MMOEAC_EnvironmentalSelection(Union,N,delta)
% Environmental selection of MMOEAC.
% Renamed from Environmental_Selection to avoid potential conflicts across algorithms.

    %% Neighborhood-based clustering method
    range = max(Union.decs,[],1) - min(Union.decs,[],1);
    r     = range * 0.1;
    C     = MMOEAC_NCM(Union.decs,r);
    K     = length(C);

    % local_A: reserve the nondominated solutions in each local cluster
    local_A = [];
    for i = 1 : K
        cluster = C(i).p;
        if length(cluster) <= delta
            continue;
        end
        [FrontNo1,~] = NDSort(Union(cluster).objs,length(cluster));
        local_A      = [local_A, cluster(FrontNo1 == 1)];
    end

    %% P: number of individuals is not less than N
    [FrontNo,~] = NDSort(Union.objs,length(Union));
    P           = find(FrontNo == 1);
    temp        = setdiff(P,local_A);
    P           = [local_A, temp];
    i = 1;
    while length(P) < N
        i    = i + 1;
        temp1 = find(FrontNo == i);
        temp2 = setdiff(temp1,local_A);
        P     = [P, temp2];
    end

    if length(P) > N
        newpop = Union(P);
        z      = min(newpop.objs,[],1);
        Z      = max(newpop.objs,[],1);
        % Normalization in the decision space
        ZZZ    = (newpop.objs - z) ./ repmat(Z - z,length(newpop),1);
        % Hierarchical clustering method
        H      = clusterdata(ZZZ,'maxclust',N,'distance','euclidean','linkage','ward');

        count = length(newpop) - N;
        for i = 1 : count
            CrowdDis = MMOEAC_Crowding(newpop.decs);
            num      = hist(H,1:N);
            % Find the most crowding clusters in the objective space
            I = find(num == max(num));
            R = [];
            for j = 1 : length(I)
                R = [R; find(H == I(j))];
            end
            % Delete the solution with minimal decision-space crowding distance
            T = find(min(CrowdDis(R)) == CrowdDis(R));
            s = randperm(length(T));
            x = R(T(s(1)));
            newpop(x) = [];
            H(x) = [];
        end
        Population = newpop;
    else
        Population = Union(P);
    end
end