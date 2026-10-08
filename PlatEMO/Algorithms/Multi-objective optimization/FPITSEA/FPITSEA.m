classdef FPITSEA < ALGORITHM
% <2023> <multi> <real> <multimodal>
% Two-stage evolutionary algorithm with fuzzy preference indicator
% p  --- 0.25 --- Ratio of FEs consumed in the first stage
% fi --- 0.1  --- Fuzzy preference threshold

%------------------------------- Reference --------------------------------
% Two-stage evolutionary algorithm with fuzzy preference indicator for
% multimodal multi-objective optimization, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [p,fi] = Algorithm.ParameterSet(0.25,0.1);

            %% Publish the active Problem so per-niche helpers can evaluate
            %% offspring against the same objective functions (FE accounting
            %% stays consistent with the main Problem).
            global PROBLEM_CGE %#ok<GLobal>
            PROBLEM_CGE = Problem;

            N          = Problem.N;
            Population = Problem.Initialization();
            [~,CD]     = FPITSEA_EnvironmentalSelectionS1(Population,N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE < p*Problem.maxFE
                    % First stage: standard fuzzy preference indicator
                    MatingPool = TournamentSelection(2,N,CD);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    Union      = [Population,Offspring];
                    [Population,CD] = FPITSEA_EnvironmentalSelectionS1(Union,N);
                else
                    % Second stage: independent evolution per niche
                    Population = FPITSEA_IndependentEvolution(Population,N,fi,Problem);
                end
            end
        end
    end
end