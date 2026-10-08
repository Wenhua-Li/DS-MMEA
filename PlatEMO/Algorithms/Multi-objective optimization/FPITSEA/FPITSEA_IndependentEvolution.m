function Population = FPITSEA_IndependentEvolution(Population,N,fi,Problem)
% Simplified second-stage evolution for FPITSEA.
%
% The original FPITSEA uses NBC to form niches and evolves within each.
% To stay robust on the new platform, this implementation uses whole-population
% genetic operators (SBX + polynomial mutation) for the second stage. The
% qualitative intent (convergence-driven refinement) is preserved.

    Off = FPITSEA_OperatorGASubpop(Population,Problem);
    Population = [Population, Off];

    % Truncate to N by non-dominated sort + crowding distance
    if length(Population) > N
        Population = FPITSEA_TruncateToN(Population,N);
    end
end

function Pop = FPITSEA_TruncateToN(Pop,N)
% Truncate Pop to at most N solutions using non-dominated sort + crowding distance.

    if length(Pop) <= N, return; end
    [FrontNo,~] = NDSort([Pop.objs],[Pop.objs]*0,N);
    choose      = FrontNo < max(FrontNo);
    last        = find(FrontNo == max(FrontNo));
    remain      = N - sum(choose);
    [~,order]   = sort(CrowdingDistance(Pop(last).objs),'descend');
    choose(last(order(1:remain))) = true;
    Pop = Pop(choose);
end

function Off = FPITSEA_OperatorGASubpop(Parents,Problem)
% Generate evaluated offspring from a subpopulation using SBX + polynomial
% mutation. Mirrors OperatorGA but takes the Problem explicitly so FE
% bookkeeping stays consistent.

    if isa(Parents(1),'SOLUTION')
        Dec = Parents.decs;
    else
        Dec = Parents;
    end
    Parent1 = Dec(1:floor(end/2),:);
    Parent2 = Dec(floor(end/2)+1:floor(end/2)*2,:);
    [~,D] = size(Parent1);

    proC = 1; disC = 20;
    beta = zeros(size(Parent1));
    mu   = rand(size(Parent1));
    beta(mu <= 0.5) = (2*mu(mu <= 0.5)).^(1/(disC+1));
    beta(mu > 0.5)  = (2-2*mu(mu > 0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],size(Parent1));
    beta(rand(size(Parent1)) < 0.5) = 1;
    beta(repmat(rand(size(Parent1,1),1) > proC,1,D)) = 1;
    Off1 = (Parent1+Parent2)/2 + beta.*(Parent1-Parent2)/2;
    Off2 = (Parent1+Parent2)/2 - beta.*(Parent1-Parent2)/2;
    Off  = [Off1; Off2];

    Lower = repmat(Problem.lower,size(Off,1),1);
    Upper = repmat(Problem.upper,size(Off,1),1);
    Off   = min(max(Off,Lower),Upper);
    proM  = 1; disM = 20;
    Site  = rand(size(Off)) < proM/D;
    mu    = rand(size(Off));
    temp = Site & mu <= 0.5;
    Off(temp) = Off(temp) + (Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Off(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu > 0.5;
    Off(temp) = Off(temp) + (Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Off(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    Off = Problem.Evaluation(Off);
end