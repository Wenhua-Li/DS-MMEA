function cand = CoSEAPropose(Population, Div, Problem, m, cfg, progress)
%CoSEAPropose Injected candidates (reflection by default, random as control).
%   Reflection is a pure geometric operation and requires no trained model
%   (design decision 5). Any failure falls back to the pure evolutionary
%   path (empty candidate set), the algorithm never dies from it.
%   cfg.proposeMode = 'random' switches the proposal distribution to uniform
%   random points in the box (A3 equal-budget control, same novelty filter,
%   same schedule, same number of evaluated candidates).

    if isfield(cfg,'proposeMode') && strcmp(cfg.proposeMode,'random')
        cand = CoSEARandomPropose(Problem, m);
    else
        cand = CoSEAReflectPropose(Population, Div, Problem, m, cfg, progress);
    end
    cand = CoSEANoveltyFilter(cand, Population, Problem, cfg);
end
