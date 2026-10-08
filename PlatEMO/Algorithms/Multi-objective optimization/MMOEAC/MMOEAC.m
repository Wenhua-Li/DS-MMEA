classdef MMOEAC < ALGORITHM
% <2019> <multi> <real> <multimodal>
% Multimodal multi-objective evolutionary algorithm with decision-space clustering

%------------------------------- Reference --------------------------------
% W. Lin, G. G. Yen, and others. MMOEAC: A multimodal multiobjective
% evolutionary algorithm with convergence and decision space clustering.
% 2019.
%------------------------------- Copyright --------------------------------
% Copyright 2019 Wu Lin
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            delta = 5;

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Mating pool: tournament selection based on decision-space density
                CrowdDis  = MMOEAC_Crowding(Population.decs);
                MatingPool = TournamentSelection(2,Problem.N,-CrowdDis);
                % Recombination
                Offspring = OperatorGA(Problem,Population(MatingPool));
                Union     = [Population,Offspring];
                % Environmental selection
                Population = MMOEAC_EnvironmentalSelection(Union,Problem.N,delta);
            end
        end
    end
end