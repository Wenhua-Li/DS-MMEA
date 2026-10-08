function fig04_dualfamily_modes()
%FIG04_DUALFAMILY_MODES Paper Fig.4: dual-family mode-count trajectories.
%   Two panels from the frozen v4.0 data (zero new runs):
%     (a) IDMPT1_NP8 (G5 sentinel): the tight family separates the 8 true
%         modes while the wide family and the v3.1 single family fuse to 1-2
%         -- the mechanism behind the 50x mode-loss repair.
%     (b) Polygon_M3_D8 (G4, D8): BOTH families fuse (medians 1-2) -- the
%         honest high-dimensional boundary (2026-09-29 erratum caliber).
%   Lines: median over 10 paired runs; bands: 25-75 percentiles.
%   Output: ../out/fig04_dualfamily_modes.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig04_dualfamily_modes"

    S = fig_style();
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    dataRoot = fullfile(root,'results','Data');
    sents = {'IDMPT1_NP8', 'Polygon_M3_D8'};
    nTrue = [8, 4];                      % true mode counts (reference line)
    nR = 10;
    nSer = 3;                            % 1=v3.1 single, 2=wide, 3=tight

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 S.w2 6.9]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    for s = 1 : 2
        nm = sents{s};
        raw = cell(1,nSer); gMin = inf;
        for r = 1 : nR
            sA = load(fullfile(dataRoot,'MSCoSEA',  sprintf('%s__r%d.mat',nm,r)));
            sB = load(fullfile(dataRoot,'MRCoSEA31',sprintf('%s__r%d.mat',nm,r)));
            raw{1}{r} = sB.injlog.nModes(:);
            raw{2}{r} = sA.injlog.nModesK(:,1);
            raw{3}{r} = sA.injlog.nModesK(:,2);
            gMin = min([gMin, numel(raw{1}{r}), numel(raw{2}{r}), numel(raw{3}{r})]);
        end
        TR = zeros(gMin,nR,nSer);
        for k = 1 : nSer
            for r = 1 : nR
                TR(:,r,k) = raw{k}{r}(1:gMin);
            end
        end
        med = squeeze(median(TR,2));
        q25 = squeeze(prctile(TR,25,2));
        q75 = squeeze(prctile(TR,75,2));
        gx  = (1:gMin).';

        ax = nexttile(tl); hold(ax,'on');
        sty = {S.A.wide, S.A.single31, S.A.tight};   % plot order: wide, single, tight
        ser = [2, 1, 3];                             % data index of each style
        for k = 1 : 3
            S.shade(gx, q25(:,ser(k)), q75(:,ser(k)), sty{k}.col);
        end
        yl = yline(ax, nTrue(s), ':', 'Color',[0.35 0.35 0.35], ...
            'LineWidth',1.0, 'HandleVisibility','off');
        hh = gobjects(1,3);
        for k = 1 : 3
            hh(k) = plot(ax, gx, med(:,ser(k)), 'LineStyle',sty{k}.ls, ...
                'Color',sty{k}.col, 'LineWidth',sty{k}.lw);
        end
        S.ax(ax);
        xlabel(ax,'Generation');
        text(ax, 0.01*gMin, nTrue(s), sprintf('  %d true modes', nTrue(s)), ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
            'FontName',S.font, 'FontSize',S.fsNote, 'Color',[0.35 0.35 0.35]);
        if s == 1
            ylabel(ax,'Number of identified modes');
            ylim(ax,[0 14]);
        else
            ylim(ax,[0 6]);
            legend(ax, hh, {sty{1}.disp, sty{2}.disp, sty{3}.disp}, ...
                'Location','northeast','FontSize',S.fsLeg,'Box','off');
        end
        S.paneltag(ax, sprintf('(%c)', char(96+s)));
    end
    S.export(fig, 'fig04_dualfamily_modes', S.w2);
    close(fig);
end
