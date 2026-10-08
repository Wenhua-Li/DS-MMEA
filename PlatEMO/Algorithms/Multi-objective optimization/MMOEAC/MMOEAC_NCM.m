function C = MMOEAC_NCM(X,r)
% Neighborhood-based clustering method (NCM) used by MMOEAC.
%
% X - matrix of decision variables (n x dim)
% r - vector of radii for each dimension (1 x dim)
%
% C - struct array, where C(k).p contains the indices of points in cluster k

    [n,dim] = size(X);
    S       = (1:n);
    D       = zeros(n,n);
    for i = 1 : dim
        D = D + (abs(X(:,i) - repmat(X(:,i)',n,1)) <= r(i));
    end

    K = 0;
    while ~isempty(S)
        K = K + 1;
        Q = [];
        C(K).p = [];
        s = randperm(length(S));
        x = S(s(1));
        Q = [Q, x];
        C(K).p = [C(K).p, x];
        while ~isempty(Q)
            ss = randperm(length(Q));
            y  = Q(ss(1));
            B  = find(D(y,:) == dim);
            T  = setdiff(B,C(K).p);
            Q  = [Q, T];
            C(K).p = [C(K).p, T];
            Q(ss(1)) = [];
        end
        S = setdiff(S,C(K).p);
    end
end