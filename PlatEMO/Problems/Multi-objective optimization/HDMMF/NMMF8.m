classdef NMMF8 < PROBLEM
% <multi> <real> <multimodal>
% Higher-dimensional multimodal multiobjective test function (HDMMF suite).
% Ported to the current platform on 2026-09-27 (see HDMMF_load for the
% path fix); objective function kept verbatim, including the original
% stale-loop-index usage of xx(i,1) inside the vectorized Q expression
% (the reference PS/PF data was generated with this exact code; see the
% problem-health-check record).

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
    end
    methods
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 12; end
            obj.lower = [zeros(1,obj.M) + 1e-10, -6.*ones(1,obj.D-obj.M)];
            obj.upper = [1.*ones(1,obj.M) - 1e-10, 6.*ones(1,obj.D-obj.M)];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            Pop = PopDec;
            M = obj.M;
            OptX = 0.2;

            obj.THETA_ = zeros(N,1);
            for i = 1:N
                obj.THETA_(i) = 2/pi*atan(Pop(i,2)./Pop(i,1));
                if Pop(i,1) == 0
                    obj.THETA_(i) = 1;
                end
            end

            xx = PopDec(:,1:M);
            dd = 2;
            x1 = PopDec(:,M+1:M+dd);       t1 = abs(sum(x1,2))./6;
            x2 = PopDec(:,M+dd+1:M+2*dd);  t2 = abs(sum(x2,2))./3;
            x3 = PopDec(:,M+2*dd+1:M+3*dd); t3 = sum(x3,2)./2;

            a = 2.*t1.^3 + 5.*t1.*t2 + 4.*t3 - 2.*t1.^2.*t3 - 18;
            b = t1 + t2.^3 + t1.*t3.^2 - 18;
            c = 8.*t1.^2 + 2.*t2.*t3 + 2.*t2.^2 + 3.*t3.^3 - 52;
            Q = (a.*b.^2.*c + a.*b.*c.^2 + b.^2 + (t1 + t2 - t3).^2 + ...
                20.*(1 - (xx(N,1).^2 + xx(N,2).^2))).^2;   % original code uses the stale loop index here

            obj.h = 20 - 20*exp(-0.2*sqrt(sum((PopDec(:,M+3*dd+1:end) - OptX).^2,2))./(D-M-3*dd)) + exp(1) ...
                - exp(sum(cos(2*pi.*(PopDec(:,M+3*dd+1:end) - OptX)),2)./(D-M-3*dd)); % P1

            T_ = zeros(N,1);
            G_ = zeros(N,M);
            for i = 1:N
                T_(i) = (1 - Pop(i,1)^2 - Pop(i,2)^2).^2 + Q(i) + obj.h(i);
                G_(i,1:M) = 1 - [ones(1,1) cumprod(sin(pi/2*obj.THETA_(i)),2)].*[cos(pi/2*obj.THETA_(i)) ones(1,1)];
            end
            PopObj = G_.*repmat((1 + T_),1,M);
        end
        function R = GetOptimum(obj,N) %#ok<INUSD>
            S = HDMMF_load('NMMF8');
            obj.POS = S.PS;
            R = S.PF;
        end
        function R = GetPF(obj) %#ok<INUSD>
            S = HDMMF_load('NMMF8');
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
