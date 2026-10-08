function fig06_rankmap()
%FIG06_RANKMAP Paper Fig.6: mean-rank heatmap of the main comparison.
%   Rows = the 8 main-comparison algorithms sorted by overall rank
%   (MSCoSEA on the first row), columns = groups G1-G6 + Overall. Each cell
%   is the Friedman average rank of the IGDX median over the 10 paired runs
%   (71 problems; identical methodology to the frozen 71-problem report and
%   to export_tab2_main.m: ties share the average rank, DNF cells take the
%   worst value). Sequential white-to-blue colormap, rank 1 = darkest; the
%   best value of every column is bold; the MSCoSEA row label is bold
%   vermillion (fixed algorithm color).
%   Output: ../out/fig06_rankmap.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig06_rankmap"

    S = fig_style();
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);       % repository root
    dataRoot = fullfile(root, 'results','Data');

    probList = rq3_problems();
    probs  = {probList.name};
    groups = {probList.group};
    nP = numel(probs);
    runs = 1:10;

    cfg = {'MSCoSEA','MSCoSEA_nomani','MRCoSEA31', ...
           'HREA','MMOEADC','CMMO','TriMOEATAR','HHCMMEA','HDMMODE','FPITSEA'};
    mainIdx = [1 4 5 6 7 8 9 10];
    dispName = {'MSCoSEA','HREA','MMOEADC','CMMO','TriMOEA-TA&R', ...
                'HHC-MMEA','HDMMODE','FPITSEA'};
    nM = numel(mainIdx);

    %% load paired runs (IGDX medians)
    MED = nan(nM, nP);
    for k = 1 : nM
        for ip = 1 : nP
            vx = [];
            for r = runs
                f = fullfile(dataRoot, cfg{mainIdx(k)}, ...
                             sprintf('%s__r%d.mat', probs{ip}, r));
                if exist(f, 'file') ~= 2, continue; end
                s = load(f, 'IGDX');
                vx(end+1) = s.IGDX;                        %#ok<AGROW>
            end
            if numel(vx) == 0, continue; end
            MED(k,ip) = median(vx);
        end
    end

    R = localRanks(MED);
    ovm = mean(R, 2);
    gTags = {'G1','G2','G3','G4','G5','G6'};
    gRm = zeros(nM, numel(gTags));
    for gi = 1 : numel(gTags)
        gRm(:,gi) = mean(R(:, strcmp(groups, gTags{gi})), 2);
    end
    V = [gRm ovm];                       % rows = algorithms, cols = groups+all
    colHeads = [{'MMF','IDMP','Local-front','Polygon','Mode-count','HDMMF'}, ...
                {'Overall'}];            % descriptive headers, no group codes
    fprintf('fig06 matrix (rows x [G1..G6 Overall], console cross-check):\n');
    disp(array2table(V, 'RowNames', dispName, 'VariableNames', ...
        [gTags, {'Overall'}]));

    %% order rows by overall rank (MSCoSEA stays row 1)
    [~, ord] = sort(ovm);
    V = V(ord, :);
    names = dispName(ord);
    isMS = strcmp(names, 'MSCoSEA');
    nRw = size(V, 1); nCl = size(V, 2);

    %% sequential white-to-blue colormap, best rank = darkest
    vLo = min(V(:)); vHi = max(V(:));
    nC = 256;
    t = linspace(0, 1, nC).';
    dark = S.C.blue * 0.60;              % deep end of the scale
    cmap = (1 - t) .* dark + t * [1 1 1];

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 18.1 7.6]);
    ax = axes('Parent',fig); hold(ax,'on');
    imagesc(ax, V, [vLo vHi]);
    set(ax,'YDir','reverse');            % hold() blocks the auto-flip
    colormap(ax, cmap);

    % thin white cell separators
    for k = 1 : nCl-1
        plot(ax, [k k]+0.5, [0.5 nRw+0.5], '-', 'Color','w', 'LineWidth',0.6);
    end
    for k = 1 : nRw-1
        plot(ax, [0.5 nCl+0.5], [k k]+0.5, '-', 'Color','w', 'LineWidth',0.6);
    end

    % cell text: column-best in bold; white text on dark cells
    for rr = 1 : nRw
        for cc = 1 : nCl
            v = V(rr,cc);
            tv = (v - vLo) / max(vHi - vLo, eps);
            cellCol = (1 - tv) .* dark + tv * [1 1 1];
            lum = 0.299*cellCol(1) + 0.587*cellCol(2) + 0.114*cellCol(3);
            if lum < 0.45, tc = [1 1 1]; else, tc = [0.12 0.12 0.12]; end
            fw = 'normal';
            if v <= min(V(:,cc)) + 1e-9, fw = 'bold'; end
            text(ax, cc, rr, sprintf('%.2f', v), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontName',S.font,'FontSize',S.fsNote,'Color',tc, ...
                'FontWeight',fw);
        end
    end

    set(ax,'XTick',1:nCl,'XTickLabel',colHeads, ...
           'YTick',1:nRw,'YTickLabel',repmat({''},nRw,1), ...
           'TickLength',[0 0],'Box','off','XAxisLocation','top');
    set(ax,'FontName',S.font,'FontSize',S.fsAx,'LineWidth',0.75);

    % row labels as normalized text objects (per-label color and weight)
    for rr = 1 : nRw
        if isMS(rr), tc = S.C.verm; fw = 'bold'; else, tc = [0.12 0.12 0.12]; fw = 'normal'; end
        yy = 1 - (rr - 0.5) / nRw;       % data row rr in normalized height
        text(ax, -0.015, yy, names{rr}, 'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','middle', ...
            'FontName',S.font,'FontSize',S.fsAx,'Color',tc,'FontWeight',fw);
    end
    xlim(ax,[0.5, nCl+0.5]);
    ylim(ax,[0.5, nRw+0.5]);            % imagesc already reversed YDir

    cb = colorbar(ax);
    cb.FontName = S.font; cb.FontSize = S.fsAx;
    cb.Label.String = 'Mean rank of IGDX median';
    cb.Label.FontName = S.font; cb.Label.FontSize = S.fsAx;
    cb.Ticks = 2:2:8;

    S.export(fig, 'fig08_rankmap', 18.1);
    close(fig);
end

function R = localRanks(M)
%localRanks Friedman-style average ranks per column (small = good).
    R = zeros(size(M));
    for ip = 1 : size(M, 2)
        tbl = M(:, ip);
        tbl(~isfinite(tbl)) = inf;
        [~, rk] = sort(tbl);
        [~, rk] = sort(rk);
        for vv = unique(tbl(:)).'
            idx = tbl == vv;
            if sum(idx) > 1
                rk(idx) = mean(rk(idx));
            end
        end
        R(:, ip) = rk;
    end
end
