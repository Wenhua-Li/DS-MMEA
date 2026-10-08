function fig07_coverage_evolution()
%fig07_coverage_evolution Paper figure: delivered coverage across the budget.
%   Covered truth modes / nModes of the delivered population at the five 20%
%   budget checkpoints, MSCoSEA versus MSCoSEA-1F, group medians. Modes are
%   found early everywhere; the differences appear in what is kept: on the
%   Polygon group both algorithms lose coverage mid-run (chaotic trajectories)
%   and MSCoSEA keeps more; on the mode-count group the single-family
%   counterpart drops a mode at the end. On the HDMMF group the 0.05
%   box-normalized membership threshold degenerates at this dimensionality,
%   so that panel shows the median distance to the true set instead (log scale).
%   Data basis: checkpoint populations of the frozen v4.2 runs (requires
%   exp/behavior/behavior_stats.m run once; reads behavior.json).
%   Output: ../out/fig07_coverage_evolution.{pdf,png}

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    J = jsondecode(fileread(fullfile(root,'results','behavior','behavior.json')));
    S = fig_style();

    gDisp = {'MMF','IDMP','Local-front','Polygon','Mode-count','HDMMF'};
    gTag  = {'G1','G2','G3','G4','G5','G6'};

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 18.1 6.4]);
    tl = tiledlayout(1,6,'TileSpacing','compact','Padding','compact');
    for gi = 1:6
        ax = nexttile(tl); hold(ax,'on');
        sel = strcmp({J.prob.group},gTag{gi});
        if gi < 6
            vM = median(J.covM(sel,:) ./ J.truthNm(sel),1);
            vF = median(J.covF(sel,:) ./ J.truthNm(sel),1);
            hM = plot(ax,0:4,vM,'-o','Color',S.C.verm,'LineWidth',1.5,'MarkerSize',4);
            hF = plot(ax,0:4,vF,'--s','Color',S.C.gray,'LineWidth',1.3,'MarkerSize',4);
            ylim(ax,[0 1.05]);
            ttl = gDisp{gi};
        else
            dM = median(J.sprM(sel,:),1);
            dF = median(J.sprF(sel,:),1);
            hM = plot(ax,0:4,dM,'-o','Color',S.C.verm,'LineWidth',1.5,'MarkerSize',4);
            hF = plot(ax,0:4,dF,'--s','Color',S.C.gray,'LineWidth',1.3,'MarkerSize',4);
            set(ax,'YScale','log');
            ttl = [gDisp{gi} ' (dist.)'];
        end
        set(ax,'FontSize',8,'Color','w','XTick',0:4, ...
            'XTickLabel',{'20','40','60','80','100'});
        title(ax,ttl,'FontSize',9);
        if gi == 1
            ylabel(ax,'covered modes / true modes','FontSize',8);
        end
    end
    set(tl,'XLabel',xlabel(tl,'budget spent (%)','FontSize',8));
    title(tl,'Delivered-set coverage at the 20% budget checkpoints (group medians)','FontSize',9);
    exportgraphics(fig,fullfile(here,'out','fig07_coverage_evolution.pdf'), ...
        'ContentType','vector');
    exportgraphics(fig,fullfile(here,'out','fig07_coverage_evolution.png'), ...
        'Resolution',300);
    fprintf('fig_style: wrote fig07_coverage_evolution.{pdf,png}\n');
end
