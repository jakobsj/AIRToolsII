%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR Tools II
% Version 2.0 of AIR Tools
% 2017-03-08
% 
% Per Christian Hansen and Jakob Sauer Jorgensen
% Dept. Applied Mathematics and Computer Science (DTU Compute)
% Technical University of Denmark
%
% Contact: pcha@dtu.dk
% 
% COPYRIGHT NOTICE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Info:
%   CONTENTS.m          - This file with overview of available files.
%   CHANGES.txt         - List of changes made in code since last version.
%
% Iterative ART Methods:
%   art                 - General interface for all Kaczmarz/ART methods.
%   kaczmarz            - Kaczmarz's method (often referred to as ART).
%   randkaczmarz        - Randomized Kaczmarz method.
%   symkaczmarz         - Symmetric Kaczmarz method.
%
% Iterative CART Methods:
%   cart                - General interface for all column-action methods.
%   columnaction        - Column-action method with cyclic column sweeps.
%
% Iterative SIRT Methods:
%   sirt                - General interface for all SIRT methods.
%   cav                 - Component Averaging (CAV) method.
%   cimmino             - Cimmino's method.
%   drop                - Diagonally Relaxed Orthogonal Projections (DROP)
%                         method.
%   landweber           - The classical Landweber method.
%   sart                - Simultaneous Algebraic Reconstruction Technique
%                         (SART) method.
%
% Training Routines:
%   train_dpme          - Training method for the stopping rules DP and ME.
%   train_relaxpar_art  - Training to determine optimal relaxation parameter
%                         for ART and CART methods.
%   train_relaxpar_sirt - Training to determine optimal relaxation parameter
%                         for SIRT method.
%
% Test Problems:
%   paralleltomo        - 2D tomography test problem, parallel beam.
%   fancurvedtomo       - 2D tomography test problem, fan beam, arc detector.
%   fanlineartomo       - 2D tomography test problem, fan beam, linear detector.
%   phantomgallery      - A collection of different 2D phantoms.
%   purge_rows          - Remove empty or very sparse rows of A and the
%                         corresponding elements of b.
%   seismictomo         - 2D seismic travel-time tomography test problem.
%   seismicwavetomo     - Similar to seismictomo but without a ray assumption.
%
% Demonstration Scripts:
%   demo_art            - Demonstrates the use of, and the results from,
%                         the ART methods.
%   demo_cart           - Demonstrates the use of, and the results from,
%                         the CART method.
%   demo_constraints    - Demonstrates the use of nonnegativity and other
%                         constraints.
%   demo_custom         - Demonstrates how to specify custom ART and SIRT
%                         methods.
%   demo_matrixfree     - Demonstrates the matrix-free mode for test problems
%                         and reconstruction methods.
%   demo_sirt           - Demonstrates the use of, and the results from,
%                         the SIRT methods.
%   demo_training       - Demonstrates the use of the training methods.
%
% Auxiliary Routines:
%   afun_matrix         - Wrap a matrix into a "matrix-free" function handle.
%   calc_relaxpar_art   - Computes the relaxation parameter for ART methods.
%   calc_relaxpar_cart  - Computes the relaxation parameter for CART methods.
%   calc_relaxpar_sirt  - Computes the relaxation parameter for SIRT methods.
%   check_inputs        - Checks inputs and sets default values.
%   check_stoprules     - Checks if stopping rule criteria are met.
%   get_mfun_dfun       - Computes matrices M and D for SIRT methods.