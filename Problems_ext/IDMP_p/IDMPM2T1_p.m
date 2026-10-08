classdef IDMPM2T1_p < PROBLEM
% <2019> <multi> <real> <multimodal>
% Parameterized imbalanced distance minimization problem (2 objectives).
%   NP modal tracks at x2 = centers(j), j = 1..NP, with centers uniformly
%   spaced in [-0.5,0.5] (centrally symmetric); on track j the Pareto set is
%   the segment x1 in [centers(j)-psize, centers(j)+psize]. Track weights
%   alternate between 1 (odd j) and a (even j): imbalanced search difficulty
%   (NOT optimality), consistent with the T2 power pattern [2, 2-a].
%   Parameters ('parameter',{a,NP,psize,localPen}): a = 3, NP = 2,
%   psize = 0.10, localPen = 0 reproduce the ORIGINAL IDMPM2T1 (all tracks
%   map to the same PF, verified against TruePSPFdata/IDMPM2T1.mat);
%   localPen = 0.01 reproduces IDMPM2T1_e bit-identically (local fronts).
%   NOTE (reconstruction 2026-09): the NP>2 generalization (linspace centers,
%   escalating weights) is a reconstruction choice - the original
%   CODES-thesis IDMP_p was lost and no NP=4/8 reference data survives.

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, G. G. Yen, Y. Nojima and N. Masuyama, Handling
% Imbalance Between Convergence and Diversity in the Decision Space in
% Evolutionary Multimodal Multi-Objective Optimization, IEEE TEVC, 2019.
%--------------------------------------------------------------------------

    properties
        a;        % imbalance base
        NP;       % number of modal tracks (equivalent Pareto sets)
        psize;    % half length of each Pareto set segment
        localPen; % local-front penalty per objective on x2>0 (0 = original)
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.a, obj.NP, obj.psize, obj.localPen] = obj.ParameterSet(3, 2, 0.10, 0);
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [-1,-1];
            obj.upper    = [1,1];
            obj.encoding = ones(1,obj.D);
        end
        %% Centers of the modal tracks (centrally symmetric)
        function c = centers(obj)
            if obj.NP == 1
                c = 0;
            else
                c = linspace(-0.5, 0.5, obj.NP);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,~]  = size(X);
            c      = obj.centers();
            w      = ones(1,obj.NP);
            w(2:2:end) = obj.a;         % alternating imbalance [1,a,1,a,...]
            Points = [c - obj.psize; c + obj.psize];
            PopObj = NaN(N,obj.M);
            for i = 1 : obj.M
                temp = abs(repmat(X(:,1),1,obj.NP) - repmat(Points(i,:),N,1));
                for j = 1 : obj.NP
                    temp(:,j) = temp(:,j) + w(j) * abs(X(:,2) - c(j));
                end
                PopObj(:,i) = min(temp,[],2);
            end
            % Local front penalty (per objective)
            if obj.localPen > 0
                index = X(:,2) > 0;
                PopObj(index,:) = PopObj(index,:) + obj.localPen;
            end
        end
        %% Generate Pareto optimal solutions (all modal tracks, per-mode cells)
        function R = GetOptimum(obj,N)
            c    = obj.centers();
            nPer = max(2, ceil(N/obj.NP));
            x1   = linspace(-obj.psize, obj.psize, nPer);
            R    = [];
            for j = 1 : obj.NP
                R = [R; c(j)+x1.', repmat(c(j),nPer,1)]; %#ok<AGROW>
            end
        end
        %% True Pareto sets per mode (for nicheRecovery)
        function modes = GetPSModes(obj,N)
            c    = obj.centers();
            nPer = max(2, ceil(N/obj.NP));
            x1   = linspace(-obj.psize, obj.psize, nPer);
            modes = cell(1,obj.NP);
            for j = 1 : obj.NP
                modes{j} = [c(j)+x1.', repmat(c(j),nPer,1)];
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj,~)
            temp = 0:0.002:2*obj.psize;
            R = [temp', 2*obj.psize-temp'];
        end
    end
end
