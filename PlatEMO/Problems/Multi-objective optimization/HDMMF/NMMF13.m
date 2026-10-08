classdef NMMF13 < PROBLEM
% <multi> <real> <multimodal>
% Higher-dimensional multimodal multiobjective test function (HDMMF suite,
% three objectives). Ported to the current platform on 2026-09-27 (see
% HDMMF_load for the path fix); objective function unchanged. GetPF keeps
% the original analytic sphere-grid form for M = 3.

%------------------------------- Reference --------------------------------
% HDMMF benchmark of: Multiobjective differential evolution for
% higher-dimensional multimodal multiobjective optimization (PDF bundled
% with the suite).
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    properties
        POS;    % Pareto optimal set (decision space), for IGDX
        Sx
        h
        THETA_
        dd
        num_of_peaks
        K
        l
    end
    methods
        function Setting(obj)
            obj.M = 3;
            obj.dd = 50;
            obj.K = 3;
            if isempty(obj.D); obj.D = 60; end
            obj.l = [];                                 % unused in the original
            obj.lower = [zeros(1,obj.K) + 1e-10, -10.*ones(1,obj.dd), ...
                         zeros(1,obj.D-obj.K-obj.dd)];
            obj.upper = [ones(1,obj.K) - 1e-10, 10.*ones(1,obj.dd), ...
                         ones(1,obj.D-obj.K-obj.dd)];
            obj.encoding = ones(1,obj.D);
            obj.num_of_peaks = 2;
        end
        function PopObj = CalObj(obj,PopDec)
            OptX = 0.2;
            [N,~] = size(PopDec);
            M = obj.M;

            TT = cumsum(PopDec(:,1:M).^2,2,'reverse');
            obj.THETA_ = zeros(N,M-1);
            for i = 1:N
                for j = 1:M-1
                    obj.THETA_(i,j) = 2/pi*atan(TT(i,j+1)./PopDec(i,j));
                    if PopDec(i,j) == 0
                        obj.THETA_(i,j) = 1;
                    end
                end
            end

            for i = obj.K+1 : 2 : obj.K+obj.dd-1
                PopDec(:,i) = (1 + cos(pi/2.*(i/(obj.K+obj.dd-1)))).*((1.5).*PopDec(:,i) - obj.lower(i)) - PopDec(:,1).*(obj.upper(i) - obj.lower(i));
            end
            for i = obj.K+2 : 2 : obj.K+obj.dd-1
                PopDec(:,i) = (1 + (i/(obj.K+obj.dd-1))).*((1.5).*PopDec(:,i) - obj.lower(i)) - PopDec(:,2).*(obj.upper(i) - obj.lower(i));
            end

            t = PopDec(:,obj.K+1:obj.K+obj.dd-1);
            n = obj.K+2:obj.K+obj.dd-1;
            for k = 1:N
                Q(k,:) = 100.*(1 - sin((1/8).*obj.num_of_peaks*pi*(t(k,1)-2))).^2 + ...
                    50.*sum(n.*(t(k,2:end) + t(k,1:end-1)).^2) + ...
                    50.*(obj.K+obj.dd).*(2.*PopDec(k,obj.K+obj.dd) - PopDec(k,1) - PopDec(k,2)).^2;
            end

            obj.h = 20 - 20*exp(-0.2*sqrt(sum((PopDec(:,obj.K+obj.dd+1:end) - OptX).^2,2))./(obj.D-(obj.K+obj.dd))) + exp(1) ...
                - exp(sum(cos(2*pi.*(PopDec(:,obj.K+obj.dd+1:end) - OptX)),2)./(obj.D-(obj.K+obj.dd)));

            T_ = zeros(N,1);
            G_ = zeros(N,M);
            for i = 1:N
                T_(i) = (1 - TT(i,1)^2).^2 + Q(i) + obj.h(i);
                G_(i,1:M) = [ones(1,1) cumprod(sin(pi/2*obj.THETA_(i,:)),2)].*[cos(pi/2*obj.THETA_(i,:)) ones(1,1)];
            end
            PopObj = G_.*repmat((1 + T_),1,M);
        end
        function R = GetOptimum(obj,N) %#ok<INUSD>
            S = HDMMF_load('NMMF13');
            obj.POS = S.PS;
            R = S.PF;
        end
        function R = GetPF(obj)
            if obj.M == 2
                R = UniformPoint(100,obj.M);
                R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'), sin(a)*sin(a'), cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
        function score = CalMetric(obj,metName,Population)
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end
