function fig07_popshow()
%FIG07_POPSHOW Paper Fig.7: final populations on four representative problems.
%   Grid: 4 problems (rows) x 4 panels (columns):
%     [MSCoSEA decision | MMOEADC decision | MSCoSEA objective | MMOEADC obj.]
%   In every panel the reference PS/PF (the same truth used for scoring, from
%   the CoSEATruePS cache) is drawn underneath in light grey, and the final
%   population of one frozen run is drawn on top in the algorithm color.
%   Shown run = the run whose IGDX is closest to that algorithm's own median
%   IGDX over its frozen runs (runs 1-10; Polygon_M4_D4 follows the 20-run
%   convention) -- the per-algorithm "median run", never a hand-picked case.
%
%   Decision space: D = 2 -> plain (x1,x2) scatter; D > 3 -> the reference PS
%   and BOTH final populations of the row are pooled, box-normalized, and ONE
%   PCA basis (2 components) is fitted on the pool; every decision panel of
%   the row is projected with that same basis. Objective space: M = 2 ->
%   (f1,f2) scatter; M > 3 -> parallel coordinates. Non-finite objective rows
%   are dropped. Panels of the same row and the same space share identical
%   axis limits. D = 3 / M = 3 fall back to 3-D scatter (not hit by these
%   four problems, kept for parity with plot_pop20_mr71).
%
%   Data (frozen, zero new evaluations):
%     results/Data/<alg>/<prob>__r<k>.mat (final PopDec/PopObj).
%   Output: ../out/fig07_popshow.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig07_popshow"

    S = fig_style();
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    dataRoot = fullfile(root,'results','Data');

    probs  = {'MMF1'; 'IDMPT1_NP4'; 'Polygon_M4_D4'; 'IDMPM2T4_e'};
    algs   = {'MSCoSEA', 'MMOEADC'};
    algCol = [S.A.MSCoSEA.col; S.A.MMOEADC.col];
    refC   = [0.75 0.75 0.75];              % reference PS/PF grey
    colT = {'MSCoSEA', 'MMOEADC', 'MSCoSEA', 'MMOEADC'; ...
            'decision space', 'decision space', 'objective space', ...
            'objective space'};   % two-line column headers: name / space

    nR = numel(probs); nC = numel(algs);
    PopD = cell(nR,nC);                     % final PopDec of the shown run
    PopO = cell(nR,nC);                     % finite final PopObj
    rStar = nan(nR,nC);                     % shown run index
    igx = nan(nR,nC); igd = nan(nR,nC);     % IGDX / IGD of the shown run
    Tr = cell(nR,1);                        % truth struct per problem

    % ---- load frozen data, pick the per-algorithm median-IGDX run ---------
    for r = 1:nR
        nm = probs{r};
        pl = rq3_problems();
        pe = pl(strcmp({pl.name},nm));
        if isempty(pe.param)
            Pro = pe.fcn('N',pe.N,'maxFE',pe.maxFE);
        else
            Pro = pe.fcn('N',pe.N,'maxFE',pe.maxFE,'parameter',pe.param);
        end
        Tr{r} = CoSEATruePS(Pro,nm);
        if strcmp(nm,'Polygon_M4_D4'), runs = 1:20; else, runs = 1:10; end

        for a = 1:nC
            vx = nan(1,runs(end));
            for k = runs
                f = fullfile(dataRoot,algs{a},sprintf('%s__r%d.mat',nm,k));
                if exist(f,'file') ~= 2, continue; end
                try
                    s1 = load(f,'IGDX'); vx(k) = s1.IGDX;
                catch
                end
            end
            cand = find(isfinite(vx));
            assert(~isempty(cand),'fig07: no completed run for %s / %s', ...
                   nm,algs{a});
            [~,jj] = min(abs(vx(cand) - median(vx(cand))));
            rStar(r,a) = cand(jj);          % ties resolved to smallest index
            f = fullfile(dataRoot,algs{a}, ...
                         sprintf('%s__r%d.mat',nm,rStar(r,a)));
            st = load(f,'PopDec','PopObj','IGDX','IGD');
            PopD{r,a} = st.PopDec;
            ob = st.PopObj;
            PopO{r,a} = ob(all(isfinite(ob),2),:);
            igx(r,a) = st.IGDX; igd(r,a) = st.IGD;
        end
    end

    % ---- one 4x4 grid ------------------------------------------------------
    fig = figure('Color','w','Visible','off','Units','centimeters', ...
                 'Position',[2 2 S.w2 19.0]);
    tl = tiledlayout(fig,nR,2*nC,'TileSpacing','compact','Padding','compact');

    axs = gobjects(nR,2*nC);
    for r = 1:nR
        T  = Tr{r};
        lb = T.lower(:).';  ub = T.upper(:).';
        sp = max(ub - lb,eps);
        usePCA = (T.D > 3);

        % decision-space shared frame of this row
        if usePCA
            refN = localSub(T.PS,800);
            pool = [localBox(refN,lb,sp); localBox(vertcat(PopD{r,:}),lb,sp)];
            mu = mean(pool,1);
            [V,E] = eig(cov(pool));
            [e,o] = sort(diag(E),'descend');
            V2  = V(:,o(1:2));
            pct = 100*e(1:2)/max(sum(diag(E)),eps);
            axD = localPad((pool - mu)*V2,0.05);
            prj = @(X)(localBox(X,lb,sp) - mu)*V2;
        else
            refN = localSub(T.PS,1000);
            axD  = localPad([lb; ub],0.05);
        end

        % objective-space shared frame of this row
        pfF  = T.PF(all(isfinite(T.PF),2),:);
        M    = size(pfF,2);
        allO = vertcat(PopO{r,:});
        axO  = localPad([min(allO,[],1); max(allO,[],1)],0.05);
        pfS  = localSub(pfF,1000);

        for a = 1:nC
            % decision panel (columns 1-2)
            ax = nexttile((r-1)*4+a); axs(r,a) = ax; hold(ax,'on');
            if usePCA
                sR = prj(refN);
                sP = prj(PopD{r,a});
                scatter(ax,sR(:,1),sR(:,2),6.5,refC,'filled', ...
                        'MarkerFaceAlpha',1);
                scatter(ax,sP(:,1),sP(:,2),7.5,algCol(a,:),'filled', ...
                        'MarkerFaceAlpha',0.45);
                xlim(ax,axD(:,1)); ylim(ax,axD(:,2));
                xlabel(ax,sprintf('PC1 (%.0f%%)',pct(1)));
                ylabel(ax,sprintf('PC2 (%.0f%%)',pct(2)));
            elseif T.D == 2
                scatter(ax,refN(:,1),refN(:,2),6.5,refC,'filled', ...
                        'MarkerFaceAlpha',1);
                scatter(ax,PopD{r,a}(:,1),PopD{r,a}(:,2),7.5, ...
                        algCol(a,:),'filled','MarkerFaceAlpha',0.45);
                xlim(ax,axD(:,1)); ylim(ax,axD(:,2));
                xlabel(ax,'x_1'); ylabel(ax,'x_2');
            else % D == 3, parity fallback
                scatter3(ax,refN(:,1),refN(:,2),refN(:,3),5,refC,'filled');
                scatter3(ax,PopD{r,a}(:,1),PopD{r,a}(:,2),PopD{r,a}(:,3), ...
                         7.5,algCol(a,:),'filled');
                xlim(ax,axD(:,1)); ylim(ax,axD(:,2)); zlim(ax,axD(:,3));
                view(ax,3);
                xlabel(ax,'x_1'); ylabel(ax,'x_2'); zlabel(ax,'x_3');
            end
            S.ax(ax);
            if r == 1
                ht = title(ax,sprintf('%s\n%s',colT{1,a},colT{2,a}));
                set(ht,'FontName',S.font,'FontSize',S.fsAx,'FontWeight', ...
                       'bold');
            end

            % objective panel (columns 3-4)
            ax = nexttile((r-1)*4+2+a); axs(r,2+a) = ax; hold(ax,'on');
            if M == 2
                scatter(ax,pfS(:,1),pfS(:,2),6.5,refC,'filled', ...
                        'MarkerFaceAlpha',1);
                scatter(ax,PopO{r,a}(:,1),PopO{r,a}(:,2),7.5, ...
                        algCol(a,:),'filled','MarkerFaceAlpha',0.45);
                xlim(ax,axO(:,1)); ylim(ax,axO(:,2));
                xlabel(ax,'f_1'); ylabel(ax,'f_2');
            elseif M == 3
                scatter3(ax,pfS(:,1),pfS(:,2),pfS(:,3),5,refC,'filled');
                scatter3(ax,PopO{r,a}(:,1),PopO{r,a}(:,2),PopO{r,a}(:,3), ...
                         7.5,algCol(a,:),'filled');
                xlim(ax,axO(:,1)); ylim(ax,axO(:,2)); zlim(ax,axO(:,3));
                view(ax,3);
                xlabel(ax,'f_1'); ylabel(ax,'f_2'); zlabel(ax,'f_3');
            else % parallel coordinates
                pfB = localSub(pfF,300);
                plot(ax,repmat((1:M)',1,size(pfB,1)),pfB', ...
                     'Color',0.85*[1 1 1],'LineWidth',0.3);
                cb = 0.45*algCol(a,:) + 0.55*[1 1 1];
                plot(ax,repmat((1:M)',1,size(PopO{r,a},1)),PopO{r,a}', ...
                     'Color',cb,'LineWidth',0.5);
                xlim(ax,[0.8, M+0.2]);
                ylim(ax,[min(axO(1,:)), max(axO(2,:))]);
                xticks(ax,1:M);
                xlabel(ax,'objective index');
                ylabel(ax,'objective value');
            end
            S.ax(ax);
            if r == 1
                ht = title(ax,sprintf('%s\n%s',colT{1,2+a},colT{2,2+a}));
                set(ht,'FontName',S.font,'FontSize',S.fsAx,'FontWeight', ...
                       'bold');
            end
        end
    end

    % ---- reserve real margins, then place the row labels ------------------
    % (tl.Position must be set AFTER all tiles exist, otherwise adding tiles
    %  re-triggers the automatic layout and silently discards it)
    tl.Position = [0.078 0.020 0.912 0.945];
    drawnow;
    axL = axes(fig,'Position',[0 0 1 1],'Visible','off');
    for r = 1:nR
        yl = axs(r,1).YLabel;
        oldU = get(yl,'Units'); set(yl,'Units','normalized');
        ext = get(yl,'Extent');             % tight bbox in axes-normalized units
        set(yl,'Units',oldU);
        ap   = get(axs(r,1),'Position');    % axes box in figure-normalized units
        xFig = ap(1) + ext(1)*ap(3);
        yFig = ap(2) + (ext(2) + ext(4)/2)*ap(4);
        xC = max(xFig - 0.021, 0.008);      % fixed gap left of the ylabel
        text(axL,xC,yFig,probs{r},'Units','normalized','Rotation',90, ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontName',S.font,'FontSize',S.fsAx,'FontWeight','bold', ...
             'Interpreter','none');
    end

    % ---- console summary (chosen runs, for the record) --------------------
    for r = 1:nR
        for a = 1:nC
            fprintf('fig07: %-14s %-8s run %2d  IGDX %.4g  IGD %.4g  |pop| %d\n', ...
                    probs{r},algs{a},rStar(r,a),igx(r,a),igd(r,a), ...
                    size(PopD{r,a},1));
        end
    end

    S.export(fig,'fig09_popshow',S.w2);
    close(fig);
end

function Y = localBox(X,lb,sp)
%localBox Normalize columns of X into the unit box.
    Y = (X - lb) ./ sp;
end

function Y = localSub(X,n)
%localSub Subsample rows of X to at most n (deterministic even spacing).
    if size(X,1) > n
        X = X(unique(round(linspace(1,size(X,1),n))),:);
    end
    Y = X;
end

function lim = localPad(X,p)
%localPad Axis limits [lo; hi] of the row range of X, padded by fraction p.
    lo  = min(X,[],1);
    hi  = max(X,[],1);
    pad = p*max(hi - lo,eps);
    lim = [lo - pad; hi + pad];
end
