function fig05_g5_modeloss()
%FIG05_G5_MODELOSS Paper Fig.5: mode-loss repair on IDMPT1_NP8 (G5 sentinel).
%   Per-run per-mode member counts of the final populations (10 paired runs x
%   8 true modes): (a) MRCoSEA (predecessor), (b) MSCoSEA. Both panels share
%   one color scale; empty cells (mode lost in that run) get an alert tint.
%   Data: frozen v4.0 runs (zero new evaluations).
%   Output: ../out/fig05_g5_modeloss.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig05_g5_modeloss"

    S = fig_style();
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    dataRoot = fullfile(root,'results','Data');
    nm  = 'IDMPT1_NP8';
    nR  = 10;
    nMd = 8;

    % true PS and per-point mode labels (IDMP block structure: equal blocks)
    rng(1);
    Pro = IDMPM2T1_p('N',200,'maxFE',10000,'parameter',{3,8,0.10,0});
    PS  = Pro.GetOptimum(20000);
    pm  = repelem(1:nMd, size(PS,1)/nMd).';
    lb = Pro.lower; ub = Pro.upper;
    PSn = (PS - lb) ./ max(ub - lb, eps);

    algs = {'MRCoSEA31','MSCoSEA'};   % data dir names; display names below
    C = zeros(nR, nMd, 2);
    for a = 1 : 2
        for r = 1 : nR
            D  = load(fullfile(dataRoot,algs{a},sprintf('%s__r%d.mat',nm,r)));
            X  = D.ckpt{end}{2};
            Xn = (X - lb) ./ max(ub - lb, eps);
            [~, j] = min(pdist2(Xn, PSn), [], 2);
            ml = pm(j);
            for k = 1 : nMd
                C(r,k,a) = sum(ml == k);
            end
        end
    end
    cMax = max(C,[],'all');
    nEmpty = squeeze(sum(C == 0, [1 2]));

    % sequential white->blue colormap; alert tint for empty cells
    nC = 128;
    seq = [linspace(1,0.00,nC).', linspace(1,0.28,nC).', linspace(1,0.55,nC).'];
    alert = [0.96 0.80 0.76];

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 13.6 6.2]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    names = {'MSCoSEA-1F','MSCoSEA'};   % panel titles (paper names)
    axs = gobjects(1,2);
    for a = 1 : 2
        ax = nexttile(tl); hold(ax,'on');
        axs(a) = ax;
        M = C(:,:,a);
        imagesc(ax, M, [0 cMax]);
        colormap(ax, seq);
        for r = 1 : nR
            for k = 1 : nMd
                v = M(r,k);
                if v == 0
                    rectangle(ax,'Position',[k-0.5 r-0.5 1 1], ...
                        'FaceColor',alert,'EdgeColor','none');
                end
                if v == 0
                    tc = [0.70 0.20 0.12]; fw = 'bold';
                elseif v/cMax > 0.55
                    tc = [1 1 1];          fw = 'normal';
                else
                    tc = [0.15 0.15 0.15]; fw = 'normal';
                end
                text(ax, k, r, sprintf('%d',v), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontName',S.font,'FontSize',S.fsNote,'Color',tc,'FontWeight',fw);
            end
        end
        xlim(ax,[0.5 nMd+0.5]); ylim(ax,[0.5 nR+0.5]);
        set(ax,'XTick',1:nMd,'YTick',1:nR,'TickLength',[0 0]);
        set(ax,'FontName',S.font,'FontSize',S.fsAx,'LineWidth',0.75,'Box','off');
        xlabel(ax,'True mode index');
        if a == 1
            ylabel(ax,'Run index (seed-paired)');
        end
        t = title(ax, sprintf('(%c) %s', char(96+a), names{a}));
        set(t,'FontName',S.font,'FontSize',S.fsAx,'FontWeight','bold');
    end

    cb = colorbar(axs(end));
    cb.Layout.Tile = 'east';
    cb.FontName = S.font; cb.FontSize = S.fsAx;
    cb.Label.String = 'Members in mode';
    cb.Label.FontName = S.font; cb.Label.FontSize = S.fsAx;

    S.export(fig, 'fig05_g5_modeloss', 13.6);
    fprintf('fig05: empty cells  (a)=%d  (b)=%d  (of %d each)\n', ...
        nEmpty(1), nEmpty(2), nR*nMd);
    close(fig);
end
