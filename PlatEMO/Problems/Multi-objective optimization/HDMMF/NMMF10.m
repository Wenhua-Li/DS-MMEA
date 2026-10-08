classdef NMMF10 < PROBLEM
% <multi> <real> <multimodal>
% Higher-dimensional multimodal multiobjective test function (HDMMF suite).
% Ported to the current platform on 2026-09-27 (see HDMMF_load for the
% path fix); objective function unchanged. GetPF keeps the original
% analytic unit-circle form.

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
            obj.M = 2;
            obj.dd = 9;
            obj.K = 3;
            if isempty(obj.D); obj.D = 16; end
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
            [N,D] = size(PopDec); %#ok<ASGLU> % D kept as in the original
            M = obj.M;

            T(:,1) = PopDec(:,3);
            T(:,2) = sqrt(PopDec(:,1).^2 + PopDec(:,2).^2);
            obj.THETA_ = zeros(N,1);
            for i = 1:N
                obj.THETA_(i) = 2/pi*atan(T(i,2)./T(i,1));
                if T(i,1) == 0
                    obj.THETA_(i) = 1;
                end
            end

            i = obj.K + 1;
            PopDec(:,obj.K+1) = (1 + cos(pi/2.*1/2)).*(PopDec(:,obj.K+1) - obj.lower(i)) ...
                - PopDec(:,1).*(obj.upper(i) - obj.lower(i));

            y = PopDec(:,obj.K+1:obj.K+obj.dd-1);
            t = cumsum(y,2);
            n = obj.K+2:obj.K+obj.dd-1;
            for k = 1:N
                Q(k,:) = 100.*(1 - sin((1/8).*obj.num_of_peaks*pi*(t(k,1)-2))).^2 + ...
                    50.*sum(n.*(2*t(k,2:end) + t(k,1:end-1)).^2) + ...
                    50.*(obj.K+obj.dd).*(2.*PopDec(k,obj.K+obj.dd) - PopDec(k,1) - PopDec(k,2)).^2;
            end

            obj.h = sum((PopDec(:,obj.K+obj.dd+1:obj.D) - OptX).^2,2);

            T_ = zeros(N,1);
            G_ = zeros(N,M);
            for i = 1:N
                T_(i) = (1 - (T(i,1)^2 + T(i,2)^2)).^2 + Q(i) + obj.h(i);
                G_(i,1:M) = [ones(1,1) cumprod(sin(pi/2*obj.THETA_(i)),2)].*[cos(pi/2*obj.THETA_(i)) ones(1,1)];
            end
            PopObj = G_.*repmat((1 + T_),1,M);
        end
        function R = GetOptimum(obj,N) %#ok<INUSD>
            S = HDMMF_load('NMMF10');
            obj.POS = S.PS;
            R = S.PF;
        end
        function R = GetPF(obj)
            if obj.M == 2
                R = UniformPoint(100,obj.M);
                R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
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
