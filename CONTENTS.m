%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR Tools
% Version 2.0
% 2017-02-28
% 
% Per Christian Hansen and Jakob Sauer Jorgensen
% Dept. Applied Mathematics and Computer Science
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
%   cart                - Columnwise version of Kaczmarz's method.
%
% Iterative SIRT Methods:
%   sirt                - General interface for all SIRT methods.
%   cav                 - Component Averaging (CAV) method.
%   cimmino             - Cimmino's method.
%   drop                - Diagonally Relaxed Orthogonal Projections (DROP) method.
%   landweber           - The classical Landweber method.
%   sart                - Simultaneous Algebraic Reconstruction Technique (SART) method.
%
% Training Routines:
%   train_dpme          - Training method for the stopping rules DP and ME.
%   train_relaxpar_art  - Training to determine optimal lambda for ART methods.
%   train_relaxpar_sirt - Training to determine optimal lambda for SIRT method.
%
% Test Problems:
%   paralleltomo        - 2D tomography test problem using a parallel beam.
%   fanbeamtomo         - 2D tomography test problem using fan beam and arc detector.
%   fanbeamtomolinear   - 2D tomography test problem using fan beam and linear detector.
%   phantomgallery      - A collection of different 2D phantoms.
%   seismictomo         - 2D seismic travel-time tomography test problem.
%   seismicwavetomo     - Similar to seismictomo but without a ray assumption.
%
% Demonstration Scripts:
%   demo_art            - Demonstrates the use of, and the results from, the ART methods.
%   demo_constraints    - Demonstrates the use of nonnegativity and other constraints.
%   demo_sirt           - Demonstrates the use of, and the results from, the SIRT methods.
%   demo_training       - Demonstrates the use of the training methods.
%
% Auxiliary Routines:
%   afun_matrix         - Wrap a matrix into a "matrix-free" function handle.
%   calc_relaxpar_art   - Computes the relaxation parameter for ART methods.
%   calc_relaxpar_cart  - Computes the relaxation parameter for CART methods.
%   calc_relaxpar_sirt  - Computes the relaxation parameter for SIRT methods.
%   check_inputs        - Checks inputs and sets default values.
%   check_stoprules     - Checks if stopping rule criteria are met.
%   get_mfun_tfun       - Computes matrices M and T for SIRT methods.
%   init_stoprules      - Initializes stopping rules.
%   rzr                 - Remove zero rows of A and the corresponding elements of b.