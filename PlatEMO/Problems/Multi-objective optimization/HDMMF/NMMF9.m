classdef NMMF9 < PROBLEM
% <multi> <real> <multimodal>
% Higher-dimensional multimodal multiobjective test function (HDMMF suite).
% Ported to the current platform on 2026-09-27 (see HDMMF_load for the
% path fix); objective function unchanged. The original file had no
% GetPF method; the platform reference-set display now uses draw_pf.

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
        h
        THETA_
        K
        q1
        q2
        Q1
        Q2
        Q
    end
    methods
        function Setting(obj)
            obj.M = 2;
            obj.K = 3;
            if isempty(obj.D); obj.D = 25; end

            obj.q1 = ceil(obj.D*(0.1));
            obj.q2 = ceil(obj.D*(0.7)) - 3;

            obj.lower = [zeros(1,obj.K + obj.q1 + obj.q2), ...
                         zeros(1,obj.D - (obj.K + obj.q1 + obj.q2))];
            obj.upper = [ones(1,obj.K), 5.*ones(1,obj.q1 + obj.q2), ...
                         ones(1,obj.D - (obj.K + obj.q1 + obj.q2))];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,PopDec)
            OptX = 0.2;
            [N,~] = size(PopDec);
            M = obj.M;
            Pop = PopDec;

            T(:,1) = PopDec(:,3);
            T(:,2) = sqrt(PopDec(:,1).^2 + PopDec(:,2).^2);
            obj.THETA_ = zeros(N,1);
            for i = 1:N
                obj.THETA_(i) = 2/pi*atan(T(i,2)./T(i,1));
                if T(i,1) == 0
                    obj.THETA_(i) = 1;
                end
            end

            for i = obj.K+obj.q1+1 : 1 : obj.K+obj.q1+obj.q2
                PopDec(:,i) = ((1 + cos(pi./2.*(i./(obj.K+obj.q1+obj.q2)))).*(Pop(:,i) - obj.lower(i)) ...
                    - obj.THETA_.*(obj.upper(i) - obj.lower(i)));
            end

            X_Q1 = PopDec(:,obj.K+1:obj.K+obj.q1);
            obj.Q1 = MMFfun(X_Q1);

            X_Q2 = PopDec(:,obj.K+obj.q1+1:obj.K+obj.q1+obj.q2);
            obj.Q2 = powell(X_Q2);

            obj.Q = obj.Q1 + obj.Q2;

            obj.h = 20 - 20*exp(-0.2*sqrt(sum((PopDec(:,obj.K+obj.q1+obj.q2+1:end) - OptX).^2,2))./(obj.D-(obj.K+obj.q1+obj.q2))) + exp(1) ...
                - exp(sum(cos(2*pi.*(PopDec(:,obj.K+obj.q1+obj.q2+1:end) - OptX)),2)./(obj.D-(obj.K+obj.q1+obj.q2)));

            T_ = zeros(N,1);
            G_ = zeros(N,M);
            for i = 1:N
                T_(i) = (1 - (T(i,1)^2 + T(i,2)^2)).^2 + obj.Q(i) + obj.h(i);
                G_(i,1:M) = [ones(1,1) cumprod(obj.THETA_(i),2)].*[1 - obj.THETA_(i) ones(1,1)];
            end
            PopObj = G_.*repmat((1 + T_),1,M);
            PopObj = real(PopObj);
        end
        function R = GetOptimum(obj,N) %#ok<INUSD>
            S = HDMMF_load('NMMF9');
            obj.POS = S.PS;
            R = S.PF;
        end
        function R = GetPF(obj) %#ok<INUSD>
            S = HDMMF_load('NMMF9');
            R = S.draw_pf;
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

function y = powell(x)
    [N1,q1] = size(x);
    y2 = zeros(N1,1);
    y1(:,1) = x(:,1);
    y1(:,2:q1+1) = x(:,1:q1);
    for i = 2:q1-1
        y2 = y2 + (y1(:,i-1) + 10.*y1(:,i)).^2 + ...
            5.*(y1(:,i+1) - y1(:,i+2)).^2 + (y1(:,i) - 2.*y1(:,i+1)).^4 + ...
            10.*(y1(:,i-1) - y1(:,i+2)).^4;
    end
    y = y2;
end

function y = MMFfun(x)
    t = sum(x(:,1:end),2);
    y1 = 2 - sin(2*pi.*(t - (1/4)));
    y = 2 - exp(-2.*log2(y1)).*(sin(pi.*t)).^6 - 1;
end
