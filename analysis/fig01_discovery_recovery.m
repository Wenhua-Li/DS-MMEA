function fig01_discovery_recovery()
%FIG01_DISCOVERY_RECOVERY Paper Fig.1: discovery-vs-maintenance motivation.
%   Mean niche-recovery rate over the evaluation budget (20%-step checkpoints)
%   on the IDMP problems (groups G2 and G5): the backbone without reflection
%   stalls, while the reflection channel recovers nearly all modes early.
%   Series: backbone only (CoSEA_off), v0.7 replica, CoSEA (v0.8 reflect),
%   MSCoSEA (v4.0, final algorithm). Frozen data only (zero new runs).
%   Output: ../out/fig01_discovery_recovery.{pdf,png}
%
%   Usage:
%     matlab -batch "addpath(genpath('PlatEMO')); addpath(genpath('Metrics_MMO')); addpath('analysis'); cd('analysis'); fig01_discovery_recovery"

    S = fig_style();
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    dataRoot = fullfile(root,'results','Data');
    probList = rq3_problems();
    cfgSel = {'CoSEA_off','CoSEA_v07','CoSEA','MSCoSEA'};
    nCkpt = 5;

    probs = probList(ismember({probList.group},{'G2','G5'}));

    % per-problem true PS (mode-partitioned), cached once
    truth = cell(numel(probs),1);
    for ip = 1 : numel(probs)
        pe = probs(ip);
        if isempty(pe.param)
            Pro = pe.fcn('N',pe.N,'maxFE',pe.maxFE);
        else
            Pro = pe.fcn('N',pe.N,'maxFE',pe.maxFE,'parameter',pe.param);
        end
        truth{ip} = CoSEATruePS(Pro,pe.name);
    end

    % acc(cfg, problem, ckpt, run) = niche recovery rate
    accR = nan(numel(cfgSel),numel(probs),nCkpt,20);
    for ic = 1 : numel(cfgSel)
        for ip = 1 : numel(probs)
            fs = dir(fullfile(dataRoot,cfgSel{ic},[probs(ip).name '__r*.mat']));
            for k = 1 : numel(fs)
                tok = regexp(fs(k).name,'__r(\d+)\.mat$','tokens','once');
                if isempty(tok), continue; end
                D = load(fullfile(fs(k).folder,fs(k).name));
                if ~isfield(D,'ckpt') || isempty(D.ckpt), continue; end
                T = truth{ip};
                for c = 1 : min(numel(D.ckpt),nCkpt)
                    accR(ic,ip,c,str2double(tok{1})) = ...
                        nicheRecovery(D.ckpt{c}{2},T.modes,T.lower,T.upper);
                end
            end
        end
    end

    fAxis = 0.2:0.2:1.0;
    groupsUsed = {'G2','G5'};
    groupDisp = {'IDMP','mode-count'};   % descriptive panel-tag names
    sty = {struct('col',S.C.gray,'ls','--','mk','s','disp','Backbone only (no reflection)'), ...
           struct('col',S.C.sky ,'ls','-','mk','^','disp','Reflection, ungated'), ...
           struct('col',S.C.blue,'ls','-','mk','o','disp','Reflection + survival-rate gating'), ...
           struct('col',S.C.verm,'ls','-','mk','d','disp','MSCoSEA (final)')};

    fig = figure('Color','w','Visible','off', ...
        'Units','centimeters','Position',[2 2 S.w2 6.9]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    for gi = 1 : 2
        sel = strcmp({probs.group},groupsUsed{gi});
        mR = squeeze(nanmean(nanmean(accR(:,sel,:,:),4),2));   % cfg x ckpt
        ax = nexttile(tl); hold(ax,'on');
        hh = gobjects(1,numel(cfgSel));
        for ic = 1 : numel(cfgSel)
            hh(ic) = plot(ax, fAxis, mR(ic,:), ...
                'LineStyle',sty{ic}.ls,'Color',sty{ic}.col, ...
                'LineWidth',S.lwData,'Marker',sty{ic}.mk, ...
                'MarkerSize',S.msData,'MarkerFaceColor',sty{ic}.col);
        end
        S.ax(ax);
        xlabel(ax,'Used budget (fraction of maxFE)');
        if gi == 1
            ylabel(ax,'Mean mode recovery rate');
        end
        ylim(ax,[0.55 1.02]); xlim(ax,[0.2 1.0]);
        S.paneltag(ax, sprintf('(%c) %s problems', char(96+gi), groupDisp{gi}));
        if gi == 2
            legend(ax, hh, cellfun(@(e)e.disp,sty,'UniformOutput',false), ...
                'Location','southeast','FontSize',S.fsLeg,'Box','off');
        end
    end
    S.export(fig, 'fig01_discovery_recovery', S.w2);
    close(fig);
end
