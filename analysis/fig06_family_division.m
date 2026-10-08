function fig06_family_division()
%fig06_family_division Paper figure: family division of labor, all groups.
%   Late-stage (last third of generations) median mode count per family and
%   per problem, pooled by benchmark group: the wide family fuses (few, large
%   modes) and the tight family separates (more, smaller modes) on every
%   group; the separation rate (share of runs whose tight-family count
%   reaches half the true mode count) delimits where the division pays.
%   Data basis: per-generation mechanism logs of the frozen v4.2 runs
%   (requires exp/behavior/behavior_stats.m to have produced
%   exp/behavior/out/behavior.json once; this script reads that bundle).
%   Output: ../out/fig06_family_division.{pdf,png}

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    J = jsondecode(fileread(fullfile(root,'results','behavior','behavior.json')));
    S = fig_style();

    gDisp = {'MMF','IDMP','Local-front','Polygon','Mode-count','HDMMF'};
    gTag  = {'G1','G2','G3','G4','G5','G6'};
    P = J.prob;

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 18.1 11.0]);
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for gi = 1:6
        sel = strcmp({P.group},gTag{gi});
        wv = [P(sel).wideLate];  tv = [P(sel).tightLate];
        sep = median([P(sel).sepRuns]);
        ax = nexttile(tl); hold(ax,'on');
        bp = boxchart(ones(1,sum(sel)),wv,'JitterOutliers','on');
        set(bp,'BoxFaceColor',S.C.blue,'BoxFaceAlpha',0.55,'BoxEdgeColor',S.C.blue,'LineWidth',0.8);
        bp = boxchart(2*ones(1,sum(sel)),tv,'JitterOutliers','on');
        set(bp,'BoxFaceColor',S.C.orang,'BoxFaceAlpha',0.55,'BoxEdgeColor',S.C.orang,'LineWidth',0.8);
        set(ax,'XTick',[1 2],'XTickLabel',{'wide','tight'},'FontSize',8);
        set(ax,'YScale','log','FontSize',8,'Color','w');
        title(ax,sprintf('%s (sep. %d%%)',gDisp{gi},round(100*sep)),'FontSize',9);
        if gi == 1 || gi == 4
            ylabel(ax,'modes identified (log scale)','FontSize',8);
        end
    end
    title(tl,'Late-stage mode count per family, by group (median over runs; box = problems, IQR; sep. = separation rate)','FontSize',9);
    exportgraphics(fig,fullfile(here,'out','fig06_family_division.pdf'), ...
        'ContentType','vector');
    exportgraphics(fig,fullfile(here,'out','fig06_family_division.png'), ...
        'Resolution',300);
    fprintf('fig_style: wrote fig06_family_division.{pdf,png}\n');
end
