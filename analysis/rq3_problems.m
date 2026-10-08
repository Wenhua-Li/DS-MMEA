function probList = rq3_problems()
%rq3_problems Canonical benchmark list of the RQ3 experiment matrix (56 problems).
%   probList = rq3_problems() returns a struct array with fields
%     name   - unique problem name (also the result file stem)
%     group  - 'G1' low-dim MMOPs | 'G2' IDMP | 'G3' MMOPLs (local PFs)
%              'G4' Polygon (high-dim, many-objective) | 'G5' mode-count sweep
%     fcn    - class handle used to construct the problem
%     param  - cell passed as 'parameter' ({} = none)
%     D, M   - decision variables, objectives (metadata, asserted against the
%              constructed problem by verify_truth.m)
%     N      - population size, maxFE - evaluation budget
%     nModes - expected number of equivalent PSs (NaN = derived from truth)
%
%   Protocol (CoMMEA, JAS 2023, Sec. 4.1; identical to the reference runner
%   scripts idmp_run.m / idmpe_run.m / mmf_run.m): N = 100*D, maxFE = 5000*D.
%   Deviation for G4 only: N = 200 fixed (as in the original polygonrun.m,
%   which uses popsize = 200 for every D), maxFE = 5000*D.
%
%   Groups G1-G4 mirror the four groups of the CoMMEA paper so that results
%   can be read side by side with its published tables; G3 excludes
%   IDMPM4T*_e because no public reference PS/PF exists for them (the
%   original idmpe_run.m also sweeps M = 2:3 only). G5 is our own mode-count
%   sweep on the parameterized IDMP (NP = 2/4/8).

    probList = struct('name',{},'group',{},'fcn',{},'param',{}, ...
        'D',{},'M',{},'N',{},'maxFE',{},'nModes',{});

    %% G1: 13 low-dimensional MMOPs (CEC2020/MMF subset, CoMMEA group 1)
    g1 = {'MMF1',2,2; 'MMF2',2,2; 'MMF3',2,2; 'MMF4',2,2; ...
          'MMF5',2,2; 'MMF6',2,2; 'MMF7',2,2; 'MMF8',2,2; ...
          'MMF9',2,2; 'MMF14',3,3; 'MMF1_e',2,2; 'MMF1_z',2,2; 'MMF14_a',3,3};
    for i = 1 : size(g1,1)
        probList = localAdd(probList,g1{i,1},'G1',str2func(g1{i,1}),{}, ...
            g1{i,2},g1{i,3},NaN);
    end

    %% G2: 12 IDMP problems (imbalanced search difficulty)
    %    Number of equivalent PSs: 2 for M = 2 and 4 for M = 3, 4. Verified
    %    two independent ways: the published reference files store the PS as
    %    nModes blocks of |PF| points (2000 = 2x1000, 8288 = 4x2072,
    %    8100 = 4x2025), and connectivity clustering of those PS samples
    %    returns exactly the same counts.
    idmpModes = [2,4,4];
    for M = 2 : 4
        for T = 1 : 4
            nm = sprintf('IDMPM%dT%d',M,T);
            probList = localAdd(probList,nm,'G2',str2func(nm),{},M,M, ...
                idmpModes(M-1));
        end
    end

    %% G3: 16 MMOPLs (8 IDMP_e + 8 MMF with local PFs)
    %    IDMP_e mode counts include the local PSs; verified the same way as
    %    above (402 = 2x201, 603 = 3x201, 1407 = 7x201, 8288 = 4x2072,
    %    16576 = 8x2072).
    idmpeModes = [2,2,3,7; 4,4,4,8];
    for M = 2 : 3
        for T = 1 : 4
            nm = sprintf('IDMPM%dT%d_e',M,T);
            probList = localAdd(probList,nm,'G3',str2func(nm),{},M,M, ...
                idmpeModes(M-1,T));
        end
    end
    g3 = {'MMF11',2,2; 'MMF12',2,2; 'MMF13',2,3; 'MMF15',3,3; ...
          'MMF15_a',3,3; 'MMF16_l1',3,3; 'MMF16_l2',3,3; 'MMF16_l3',3,3};
    for i = 1 : size(g3,1)
        probList = localAdd(probList,g3{i,1},'G3',str2func(g3{i,1}),{}, ...
            g3{i,2},g3{i,3},NaN);
    end

    %% G4: 9 Polygon problems (M objectives/vertices x D decision variables)
    %    row = col = 2 -> 4 polygons -> 4 equivalent Pareto sets (modes).
    %    N = 200 fixed (original polygonrun.m), maxFE = 5000*D.
    for M = [3,4,6]
        for D = [4,8,10]
            nm = sprintf('Polygon_M%d_D%d',M,D);
            probList(end+1) = struct('name',nm,'group','G4','fcn',@Polygon_p, ...
                'param',{{M,D,2,2}},'D',D,'M',M,'N',200, ...
                'maxFE',5000*D,'nModes',4); %#ok<AGROW>
        end
    end

    %% G5: 6 parameterized IDMP problems (mode-count sweep NP = 2/4/8)
    for NP = [2,4,8]
        nm = sprintf('IDMPT1_NP%d',NP);
        probList(end+1) = struct('name',nm,'group','G5','fcn',@IDMPM2T1_p, ...
            'param',{{3,NP,0.10,0}},'D',2,'M',2,'N',200, ...
            'maxFE',10000,'nModes',NP); %#ok<AGROW>
    end
    for NP = [2,4,8]
        nm = sprintf('IDMPT2_NP%d',NP);
        probList(end+1) = struct('name',nm,'group','G5','fcn',@IDMPM2T2_p, ...
            'param',{{0.4,NP,0.10,0}},'D',2,'M',2,'N',200, ...
            'maxFE',10000,'nModes',NP); %#ok<AGROW>
    end

    %% G6: 15 HDMMF problems (higher-dimensional MMO, NMMF1-15; 2026-09-27)
    %    Reference PS/PF: ReferenceData/Reference_PSPF_data/NMMF*_...mat
    %    (authors' data; mode split by connectivity clustering). Protocol
    %    follows the high-dimensional precedent of G4: N = 200 fixed
    %    (N = 100*D would give populations up to 6000), maxFE = 5000*D.
    g6 = {'NMMF1',7,2;  'NMMF2',10,2; 'NMMF3',10,2; 'NMMF4',10,2; ...
          'NMMF5',16,2; 'NMMF6',12,2; 'NMMF7',16,2; 'NMMF8',12,2; ...
          'NMMF9',25,2; 'NMMF10',16,2;'NMMF11',18,2;'NMMF12',18,2; ...
          'NMMF13',60,3;'NMMF14',10,2;'NMMF15',30,2};
    for i = 1 : size(g6,1)
        probList(end+1) = struct('name',g6{i,1},'group','G6', ...
            'fcn',str2func(g6{i,1}),'param',{{}}, ...
            'D',g6{i,2},'M',g6{i,3},'N',200, ...
            'maxFE',5000*g6{i,2},'nModes',NaN); %#ok<AGROW>
    end
end

function probList = localAdd(probList,name,group,fcn,param,M,D,nModes)
    probList(end+1) = struct('name',name,'group',group,'fcn',fcn, ...
        'param',{param},'D',D,'M',M,'N',100*D,'maxFE',5000*D,'nModes',nModes);
end
