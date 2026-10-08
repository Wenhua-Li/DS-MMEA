classdef IDMPM2T2_p < PROBLEM
% <2019> <multi> <real> <multimodal>
% Parameterized imbalanced distance minimization problem, steep variant.
%   Same layout as IDMPM2T1_p but with steep track penalties:
%   100*|x2-centers(j)|^pwr(j), pwr(j) = 2 (odd j) or 2-a (even j).
%   Parameters ('parameter',{a,NP,psize,localPen}): a = 0.4, NP = 2,
%   psize = 0.10, localPen = 0 reproduce the ORIGINAL IDMPM2T2 (equivalent
%   tracks); localPen = 0.01 reproduces IDMPM2T2_e bit-identically.
%   NOTE (reconstruction 2026-09): the NP>2 generalization is a
%   reconstruction choice - the original CODES-thesis IDMP_p was lost.

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, G. G. Yen, Y. Nojima and N. Masuyama, Handling
% Imbalance Between Convergence and Diversity in the Decision Space in
% Evolutionary Multimodal Multi-Objective Optimization, IEEE TEVC, 2019.
%--------------------------------------------------------------------------

    properties
        a;        % shape parameter of the even-track penalty
        NP;       % number of modal tracks
        psize;    % half length of each Pareto set segment
        localPen; % local-front penalty per objective on x2>0 (0 = original)
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.a, obj.NP, obj.psize, obj.localPen] = obj.ParameterSet(0.4, 2, 0.10, 0);
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
            pwr    = 2*ones(1,obj.NP);
            pwr(2:2:end) = 2 - obj.a;
            Points = [c - obj.psize; c + obj.psize];
            PopObj = NaN(N,obj.M);
            for i = 1 : obj.M
                temp = abs(repmat(X(:,1),1,obj.NP) - repmat(Points(i,:),N,1));
                for j = 1 : obj.NP
                    temp(:,j) = temp(:,j) + 100 * abs(X(:,2) - c(j)).^pwr(j);
                end
                PopObj(:,i) = min(temp,[],2);
            end
            % Local front penalty (per objective)
            if obj.localPen > 0
                index = X(:,2) > 0;
                PopObj(index,:) = PopObj(index,:) + obj.localPen;
            end
        end
        %% Generate Pareto optimal solutions (all modal tracks)
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
