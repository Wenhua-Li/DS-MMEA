function fig02_scale_dilemma()
%FIG02_SCALE_DILEMMA Paper Fig.2: the scale dilemma (semi-schematic).
%   Two synthetic scenarios x two clustering radii, same chaining rule:
%     row 1 "curved modes":    wide radius keeps each curve intact,
%                              tight radius shatters it into fragments;
%     row 2 "sparse colonies": wide radius fuses the two colonies,
%                              tight radius separates them.
%   Wide radius = 0.08*sqrt(d) (fixed); tight radius = min(0.08*sqrt(d),
%   2*Q90 of 1-NN distances) -- the two hypotheses MSCoSEA runs in parallel.
%   Points are colored by cluster id; singleton noise points in gray.
%   Deterministic toy data (fixed seed), no optimization runs involved.
%   Output: ../out/fig02_scale_dilemma.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig02_scale_dilemma"

    S = fig_style();
    rng(20260929);

    %% scenario point clouds in [0,1]^2
    % row 1: two curved modes (arcs) with uneven coverage (clumps + gaps),
    %        like a real evolving population on a curved Pareto set
    tc = linspace(0.22*pi, 0.78*pi, 7);
    Xa = [];
    for sgn = [1 -1]
        cx = 0.50 + sgn*0.25;
        for k = 1 : numel(tc)
            tt = tc(k) + 0.014*pi*randn(7,1);
            Xa = [Xa; cx + sgn*0.22*cos(tt), 0.52 + 0.30*sin(tt) + 0.004*randn(7,1)]; %#ok<AGROW>
        end
    end
    X1 = Xa;

    % row 2: three sparse colonies (the third is weak), gaps chosen so the
    %        wide radius chains across them while the tight radius cannot
    rng(11);   % dedicated sub-seed: wide fuses all three, tight separates
    X2 = [0.36 + 0.014*randn(13,1), 0.50 + 0.014*randn(13,1); ...
          0.50 + 0.014*randn(13,1), 0.52 + 0.014*randn(13,1); ...
          0.64 + 0.016*randn( 6,1), 0.49 + 0.016*randn( 6,1)];

    scen = {X1, X2};
    rowLab = {'Curved modes', 'Sparse colonies'};
    colLab = {'Wide radius (0.08\surd{\itD})', 'Tight radius (min(0.08\surd{\itD}, 2\timesQ_{90}))'};
    vfmt = {{'%d modes (intact)', '%d fragments (shattered)'}, ...
            {'%d cluster (fused)', '%d modes (separated)'}};

    d = 2;
    rWide = 0.08*sqrt(d);
    clusCols = [S.C.blue; S.C.verm; S.C.green; S.C.purp; S.C.orang; ...
                S.C.sky; S.C.yell; S.C.gray];

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 12.6 11.8]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for i = 1 : 2
        X = scen{i};
        rTight = min(rWide, 2*quantile(nnDist(X), 0.90));
        for j = 1 : 2
            if j == 1, r = rWide; else, r = rTight; end
            lab = chainCluster(X, r);
            ax = nexttile(tl); hold(ax,'on');
            u = unique(lab);
            big = u(arrayfun(@(k) sum(lab==k) >= 2, u));   % clusters of >=2
            ni = 0;
            for k = reshape(u,1,[])
                sel = lab == k;
                if ismember(k, big)
                    ni = ni + 1;
                    col = clusCols(mod(ni-1, size(clusCols,1)) + 1, :);
                else
                    col = [0.60 0.60 0.60];                 % noise
                end
                plot(ax, X(sel,1), X(sel,2), '.', 'Color',col, 'MarkerSize',7);
            end
            S.ax(ax); set(ax,'YGrid','off');
            xlim(ax,[0 1]); ylim(ax,[0.1 0.95]);
            set(ax,'XTick',[],'YTick',[]);
            % verdict chip at the bottom (count computed from the clustering)
            nBig = numel(big);
            txt = sprintf(vfmt{i}{j}, nBig);
            ok = (i==1 && j==1) || (i==2 && j==2);
            if ok, tc = [0.0 0.45 0.25]; else, tc = [0.75 0.20 0.10]; end
            text(ax, 0.5, 0.13, sprintf('%s',txt), ...
                'HorizontalAlignment','center','FontName',S.font, ...
                'FontSize',S.fsNote,'FontWeight','bold','Color',tc);
            if i == 1
                t = title(ax, colLab{j});
                set(t,'FontName',S.font,'FontSize',S.fsAx,'FontWeight','bold');
            end
            if j == 1
                ylabel(ax, rowLab{i});
            end
            S.paneltag(ax, sprintf('(%c)', char(96 + (i-1)*2 + j)));
        end
    end
    S.export(fig, 'fig02_scale_dilemma', S.w1);
    close(fig);
end

function d = nnDist(X)
%nnDist Distance to the nearest neighbor for each row of X.
    D = pdist2(X, X);
    D(1:size(D,1)+1:end) = inf;
    d = min(D, [], 2);
end

function lab = chainCluster(X, r)
%chainCluster Connected components of the r-neighborhood graph (chaining).
    n = size(X,1);
    D = pdist2(X, X);
    adj = D <= r;
    lab = zeros(n,1);
    cid = 0;
    for i = 1 : n
        if lab(i) ~= 0, continue; end
        cid = cid + 1;
        % BFS flood fill
        queue = i; lab(i) = cid;
        while ~isempty(queue)
            v = queue(1); queue(1) = [];
            nb = find(adj(v,:) & (lab == 0).');
            lab(nb) = cid;
            queue = [queue nb]; %#ok<AGROW>
        end
    end
end
