classdef DSMMEA < MSCoSEA
%DSMMEA Paper-named entry point of DS-MMEA.
%   DS-MMEA (Dual-Scale Multimodal Multiobjective Evolutionary Algorithm) is
%   the name used in the TEVC paper "Dual-Scale Co-evolutionary Algorithm
%   for Multimodal Multiobjective Optimization". The reference
%   implementation lives in MSCoSEA.m (development codename); this subclass
%   adds no code of its own, it only lets fresh experiments call the
%   algorithm by its paper name:
%
%       platemo('algorithm', {@DSMMEA}, 'problem', @MMF1, 'N', 200, 'maxFE', 10000);
%
%   Ablation switches are in MSCoSEAConfig.m (default = paper
%   configuration). See the "Naming" section of README.md.
end
