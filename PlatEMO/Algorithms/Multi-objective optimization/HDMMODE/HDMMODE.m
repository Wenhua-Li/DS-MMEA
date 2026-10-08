classdef HDMMODE < ALGORITHM
% <multi> <real> <multimodal>
% Multiobjective differential evolution for higher-dimensional multimodal
% multiobjective optimization (dual-subpopulation coevolution with decision
% space niche clustering).
%
% Ported from the authors' original code (old PlatEMO 1.x API) to the
% current platform on 2026-09-27. Changes are documented in the file
% header of each helper and in the dev-plan:
%   1. All helpers renamed with the HDMMODE_ prefix (genpath safety).
%   2. The bundled copy of NDSort.m was dropped; the platform NDSort has
%      the same interface and semantics.
%   3. Environmental_Selection.m / _2.m were merged logic-preserving with
%      a defensive fix for an unreachable branch (see HDMMODE_EnvSelect1).
%   4. Penalty_Choose.m had a linear-index bug (matrix find result used as
%      row index); fixed to unique row indices (HDMMODE_PenaltyChoose).
%
%------------------------------- Reference --------------------------------
% Multiobjective differential evolution for higher-dimensional multimodal
% multiobjective optimization (the HDMMODE paper, bundled with the HDMMF
% benchmark suite; full citation to be completed from the PDF).
% The dual-population/archive template originates from the PlatEMO CMO
% coevolutionary framework: Y. Tian, T. Zhang, J. Xiao, X. Zhang, and Y. Jin,
% A coevolutionary framework for constrained multi-objective optimization
% problems, IEEE Transactions on Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)

            % Fixed internal hyper-parameters (as in the original code)
            delta            = 3;          % min niche size kept by clustering
            fk               = 1;          % Minkowski order for distances
            K                = 2;          % number of subpopulations
            last_gen         = 100;        % stagnation detection window
            change_threshold = 0.01;       % stagnation detection threshold
            Flag             = 0;          % 1 after stagnation is detected
            gen              = 1;

            change_rate = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);

            Population = Problem.Initialization();
            NND = floor(Problem.N/K);

            Populations = cell(1,K);
            index = randperm(floor(Problem.N/K)*K);
            temp  = reshape(index,K,floor(Problem.N/K));
            MaxGen = Problem.maxFE./ Problem.N; %#ok<NASGU> % kept as in the original

            for i = 1 : K
                Populations{i} = Population(temp(i,:));
            end
            Archive = cell(1,K);
            rank = 1:K;

            for i = 1 : K
                [Front_number,~] = NDSort(Populations{rank(i)}.objs,inf);
                inter_Archive = [Populations{rank(i)}(Front_number==1),Archive{rank(i)}];
                [a_Front_number,~] = NDSort(inter_Archive.objs,inf);
                Archive{rank(i)} = inter_Archive(a_Front_number==1);
                clear a_Front_number Front_number

                while length(Archive{rank(i)}) > 1*length(Populations{rank(i)})
                    Del = HDMMODE_Truncation2( Archive{rank(i)}.objs, ...
                        length(Archive{rank(i)}) - 1*length(Populations{rank(i)}) );
                    Archive{rank(i)} = Archive{rank(i)}(Del');
                    clear Del
                end
            end

            while Algorithm.NotTerminated(Population)
                for i = 1 : K
                    if i == 1
                        % main subpopulation: niche-clustering selection
                        change_rate = HDMMODE_Normalization2(Population,change_rate,gen);
                        CrowdDis{rank(i)} = HDMMODE_Crowding(Populations{rank(i)}.decs,fk);
                        MatingPool = TournamentSelection(2,2*length(Populations{rank(i)}),-CrowdDis{rank(i)});
                        Offspring  = OperatorDE(Problem, Populations{rank(i)}, ...
                            Populations{rank(i)}(MatingPool(1:end/2)), ...
                            Populations{rank(i)}(MatingPool(end/2+1:end)));

                        if Flag == 0
                            if gen > last_gen && HDMMODE_Convertion2(change_rate,gen,last_gen,change_threshold)
                                Flag = 1;
                            end
                        end

                        if Flag == 1
                            % stagnation detected: halve the subpopulation
                            % budget and keep it there (as in the original)
                            if NND ~= floor(Problem.N/K)/2
                                NND = NND - 1;
                            end
                            Populations{rank(i)} = [Populations{rank(i)},Offspring];
                            Populations{rank(i)} = HDMMODE_EnvSelect1(Populations{rank(i)}, ...
                                NND,Problem.upper,Problem.lower,delta,fk);
                        else
                            Populations{rank(i)} = [Populations{rank(i)},Offspring];
                            Populations{rank(i)} = HDMMODE_EnvSelect1(Populations{rank(i)}, ...
                                floor(Problem.N/K),Problem.upper,Problem.lower,delta,fk);
                        end
                        [Front_number,~] = NDSort(Populations{rank(i)}.objs,inf);
                        inter_Archive = [Populations{rank(i)}(Front_number==1),Archive{rank(i)}];
                        [a_Front_number,~] = NDSort(inter_Archive.objs,inf);
                        Archive{rank(i)} = inter_Archive(a_Front_number==1);
                        clear a_Front_number Front_number Offspring

                        while length(Archive{rank(i)}) > 1*length(Populations{rank(i)})
                            Del = HDMMODE_Truncation2( Archive{rank(i)}.objs, ...
                                length(Archive{rank(i)}) - 1*length(Populations{rank(i)}) );
                            Archive{rank(i)} = Archive{rank(i)}(Del');
                            clear Del
                        end
                    else
                        % helper subpopulation: penalty-based selection
                        % against the main subpopulation's archive
                        CrowdDis{rank(i)} = HDMMODE_Crowding(Populations{rank(i)}.decs,fk);
                        MatingPool = TournamentSelection(2,2*length(Populations{rank(i)}),-CrowdDis{rank(i)});
                        Offspring  = OperatorDE(Problem, Populations{rank(i)}, ...
                            Populations{rank(i)}(MatingPool(1:end/2)), ...
                            Populations{rank(i)}(MatingPool(end/2+1:end)));

                        Populations{rank(i)} = [Populations{rank(i)},Offspring];

                        [All_front,~] = NDSort(Archive{rank(1:i-1)}.objs,inf);
                        Archive_for_Compare = Archive{rank(1:i-1)}(All_front==1);

                        Range = max([Populations{rank(1)}.decs;Populations{rank(K)}.decs],[],1) - ...
                            min([Populations{rank(1)}.decs;Populations{rank(K)}.decs],[],1);
                        rr = Range*0.2;
                        Penalty_index = HDMMODE_PenaltyChoose(Populations{rank(i)},Archive_for_Compare,rr);

                        Populations{rank(i)} = HDMMODE_EnvSelect2(Populations{rank(i)}, ...
                            2*floor(Problem.N/K) - NND, Problem.upper,Problem.lower, ...
                            delta,fk,Penalty_index);
                    end
                end
                Population = [Populations{:}];
                gen = gen + 1;
            end
        end
    end
end
