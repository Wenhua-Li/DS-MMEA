function S = fig_style()
%FIG_STYLE Shared journal-figure style layer for the TEVC paper (MSCoSEA v4.0).
%   S = FIG_STYLE() returns the Okabe-Ito palette, the fixed algorithm->style
%   map (same mapping in every figure of the paper and the SI), font sizes,
%   standard figure widths in centimeters, and helper function handles:
%     S.ax(hAx)                apply axis cosmetics (fonts, box, light y-grid)
%     S.paneltag(hAx,'(a)')    bold panel tag at the top-left corner
%     S.shade(x,yLo,yHi,col)   translucent band (alpha 0.2), returns patch
%     S.export(hFig,name,wCm)  vector PDF + 300 dpi PNG into ../out
%
%   Every figure script figNN_*.m must go through this layer. See STYLE.md.
%
%   Usage:
%     S = fig_style();

    %% Okabe-Ito palette (colorblind safe)
    C.verm  = [0.835 0.369 0.000];   % vermillion  #D55E00  -> MSCoSEA
    C.blue  = [0.000 0.447 0.741];   % blue        #0072B2  -> MMOEADC / wide family
    C.sky   = [0.337 0.706 0.914];   % sky blue    #56B4E9  -> HREA
    C.green = [0.000 0.620 0.451];   % bluish green#009E73  -> CMMO
    C.orang = [0.902 0.624 0.000];   % orange      #E69F00  -> TriMOEA-TA&R / tight family
    C.yell  = [0.941 0.894 0.259];   % yellow      #F0E442  -> HHC-MMEA
    C.purp  = [0.800 0.475 0.655];   % reddish purp#CC79A7  -> HDMMODE
    C.gray  = [0.50  0.50  0.50 ];   % gray                 -> FPITSEA / predecessor
    S.C = C;

    %% fixed algorithm -> (color, line, marker, display name) map
    A = struct();
    A.MSCoSEA    = mk(C.verm , '-' , 'o', 1.5, 'MSCoSEA');
    A.MRCoSEA31  = mk(C.gray , '--', 's', 1.2, 'MSCoSEA-1F (predecessor)');
    A.MMOEADC    = mk(C.blue , '-' , 's', 1.2, 'MMOEA/DC');
    A.HREA       = mk(C.sky  , '-' , '^', 1.2, 'HREA');
    A.CMMO       = mk(C.green, '-' , 'd', 1.2, 'CMMO');
    A.TriMOEATAR = mk(C.orang, '-' , 'v', 1.2, 'TriMOEA-TA&R');
    A.HHCMMEA    = mk(C.yell , '-' , '>', 1.2, 'HHC-MMEA');
    A.HDMMODE    = mk(C.purp , '-' , 'p', 1.2, 'HDMMODE');
    A.FPITSEA    = mk(C.gray , ':' , 'h', 1.2, 'FPITSEA');
    % family roles inside MSCoSEA (mechanism figures)
    A.wide       = mk(C.blue , '-' , 'none', 1.2, 'MSCoSEA wide family');
    A.tight      = mk(C.orang, '-' , 'none', 1.5, 'MSCoSEA tight family');
    A.single31   = mk(C.gray , '--', 'none', 1.2, 'MSCoSEA-1F (single family)');
    S.A = A;

    %% typography and geometry (IEEE two-column)
    S.font    = 'Arial';
    S.fsAx    = 8;      % axis labels and ticks
    S.fsLeg   = 7.5;    % legend
    S.fsNote  = 7.5;    % in-axes annotations
    S.w1      = 8.8;    % single-column width (cm)
    S.w2      = 18.1;   % double-column width (cm)
    S.lwData  = 1.2;
    S.lwMain  = 1.5;
    S.msData  = 4.5;
    S.gridCol = [0.85 0.85 0.85];

    %% helper handles (local functions below)
    S.ax       = @localAx;
    S.paneltag = @localPanelTag;
    S.shade    = @localShade;
    S.export   = @localExport;
end

function e = mk(col, ls, mkk, lw, dispName)
%mk Constructor for one algorithm style entry.
    e = struct('col',col,'ls',ls,'mk',mkk,'lw',lw,'disp',dispName);
end

function localAx(hAx)
%localAx Apply the shared axis cosmetics.
    S = fig_style_const();
    set(hAx,'FontName',S.font,'FontSize',S.fsAx,'LineWidth',0.75, ...
        'Color','w', ...
        'Box','off','TickDir','out','Layer','top', ...
        'YGrid','on','GridColor',S.gridCol,'GridAlpha',0.6, ...
        'XGrid','off');
end

function localPanelTag(hAx, tag)
%localPanelTag Bold panel tag "(a)" at the top-left corner of an axis.
    S = fig_style_const();
    text(hAx, 0.02, 0.96, tag, 'Units','normalized', ...
        'FontName',S.font,'FontSize',S.fsAx,'FontWeight','bold', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'BackgroundColor','w','EdgeColor','none','Margin',1.5);
end

function h = localShade(x, yLo, yHi, col)
%localShade Translucent percentile band (alpha 0.2).
    x = x(:); yLo = yLo(:); yHi = yHi(:);
    h = fill([x; flipud(x)], [yLo; flipud(yHi)], col, ...
        'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility','off');
end

function localExport(hFig, name, wCm)
%localExport Vector PDF + 300 dpi PNG preview into ../out of this file's dir.
    here = fileparts(mfilename('fullpath'));
    outDir = fullfile(here, 'out');
    if exist(outDir,'dir') ~= 7
        mkdir(outDir);
    end
    set(hFig,'Units','centimeters','Color','w');
    pos = get(hFig,'Position');
    set(hFig,'Position',[pos(1) pos(2) wCm pos(4)]);
    set(hFig,'PaperPositionMode','auto');
    exportgraphics(hFig, fullfile(outDir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(hFig, fullfile(outDir,[name '.png']), 'Resolution',300);
    fprintf('fig_style: wrote %s.{pdf,png} (%.1f cm wide) to %s\n', name, wCm, outDir);
end

function S = fig_style_const()
%fig_style_const Constants needed inside helpers (avoids nested handles).
    S.font = 'Arial';
    S.fsAx = 8;
    S.gridCol = [0.85 0.85 0.85];
end
