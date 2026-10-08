classdef NMMF14 < PROBLEM
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
        K
        q1
        q2
        Q1
        Q2
        Q3
        Q
        q
    end
    methods
        function Setting(obj)
            obj.M = 2;
            obj.K = 1;
            obj.q1 = 4;
            obj.q2 = 4;
            obj.q = 4;
            if isempty(obj.D); obj.D = 10; end
            obj.lower = [zeros(1,obj.K), -10.*ones(1,obj.D-obj.K)];
            obj.upper = [ones(1,obj.K),  10.*ones(1,obj.D-obj.K)];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,PopDec)
            OptX = 0.2;
            Pop = PopDec;

            for i = obj.K+1 : 2 : obj.K+obj.q1
                PopDec(:,i) = (1 + cos(pi./2.*(i./(obj.K+obj.q1)))).*(2.*Pop(:,i) - obj.lower(i)) ...
                    - Pop(:,1).*(obj.upper(i) - obj.lower(i))/2;
            end
            for i = obj.K+2 : 2 : obj.K+obj.q1
                PopDec(:,i) = (1 + i./(obj.K+obj.q1)).*(2.*Pop(:,i) - obj.lower(i)) ...
                    - Pop(:,1).*(obj.upper(i) - obj.lower(i))/2;
            end
            for i = obj.K+obj.q1+1 : 2 : obj.K+obj.q1+obj.q2
                PopDec(:,i) = (1 + i./(obj.K+obj.q1+obj.q2)).*(2.*Pop(:,i) - obj.lower(i)) ...
                    - ((Pop(:,obj.K+1) - obj.lower(obj.K+1))./(obj.upper(obj.K+1) - obj.lower(obj.K+1))).*(obj.upper(i) - obj.lower(i));
            end
            for i = obj.K+obj.q1+2 : 2 : obj.K+obj.q1+obj.q2
                PopDec(:,i) = (1 + cos(pi./2.*(i./(obj.K+obj.q1+obj.q2)))).*(2.*Pop(:,i) - obj.lower(i)) ...
                    - ((Pop(:,obj.K+1) - obj.lower(obj.K+1))./(obj.upper(obj.K+1) - obj.lower(obj.K+1))).*(obj.upper(i) - obj.lower(i));
            end

            X_Q1 = PopDec(:,obj.K+1:obj.K+obj.q1);
            obj.Q1 = Pathological(X_Q1);
            X_Q2 = PopDec(:,obj.K+obj.q1+1:obj.K+obj.q1+obj.q2);
            obj.Q2 = SineWave(X_Q2,obj);
            obj.Q3 = sum((PopDec(:,obj.K+obj.q1+obj.q2+1:obj.D) - OptX).^2,2);
            G = 1 + obj.Q1 + obj.Q2 + obj.Q3;

            h = 1 - (PopDec(:,1)./G).^2 - (PopDec(:,1)./G).*sin(2*pi.*obj.q.*PopDec(:,1));
            PopObj(:,1) = PopDec(:,1);
            PopObj(:,2) = G.*h;
        end
        function R = GetOptimum(obj,N) %#ok<INUSD>
            S = HDMMF_load('NMMF14');
            obj.POS = S.PS;
            R = S.PF;
        end
        function R = GetPF(obj) %#ok<INUSD>
            S = HDMMF_load('NMMF14');
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

function y = SineWave(x,obj)
    [N1,q1] = size(x); %#ok<NASGU> % N1 kept as in the original
    y = zeros(N1,1);
    for i = 1:q1-1
        j1 = obj.K+obj.q1+i; j2 = obj.K+obj.q1+obj.q2;
        y = y + ((abs(x(:,i)) - j1./j2).^2 + 10.*(x(:,i+1) - ((j1+1)./j1).*x(:,i)).^2).^(1/2).* ...
            (sin((50.*((abs(x(:,i)) - j1./j2).^2 + 10.*(x(:,i+1) - ((j1+1)./j1).*x(:,i)).^2).^(1/2))).^2 + 100);
    end
end

function y = Pathological(x)
    [N1,q1] = size(x); %#ok<NASGU> % N1 kept as in the original
    y = zeros(N1,1);
    r = r_cal(x(:,1));
    for i = 2:q1
        if i == q1
            x(:,i+1) = x(:,1);
        end
        y1 = 50.*sin(sqrt(100.*(x(:,i+1) - x(:,i)).^2 + (x(:,i) - x(:,i-1)).^2)).^2 + 100*r;
        y2 = 1 + 1.*((x(:,i+1) - x(:,i)).^2 - 2.*(x(:,i+1) - x(:,i)).*(x(:,i) - x(:,i-1)) + (x(:,i) - x(:,i-1)).^2).^2;
        y = y + y1./y2;
    end
end

function r = r_cal(x)
    n = size(x,1);
    for jj = 1:n
        if x(jj,1) > 5
            r(jj,:) = 25.*(x(jj,1)-6).^2;
        elseif x(jj,1) < -5
            r(jj,:) = 25./4.*(x(jj,1)+7).^2;
        else
            r(jj,:) = (24/25).*x(jj,1).^2 + 1;
        end
    end
end
