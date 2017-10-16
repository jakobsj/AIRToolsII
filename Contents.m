%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR Tools II
% Updated and expanded version of the AIR Tools package.
% October 2017
% 
% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.
% 
% Contact: pcha@dtu.dk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Info:
%   Contents            - This file with an overview of all available files
%   AIRToolsII_setup    - Function that sets up search paths to AIR Tools II
%   CHANGES.txt         - List of changes made in code since last versions
%   LICENSE.txt         - License file
%   README.txt          - README file with overview and install guide
%
% Iterative ART Methods:
%   art                 - General interface for all Kaczmarz/ART methods
%   kaczmarz            - Kaczmarz's method with cyclic row sweep (often referred to as ART)
%   randkaczmarz        - Kaczmarz's method with random row selection
%   symkaczmarz         - Kaczmarz's method with symmetric (top-bottom-top) row sweep
%
% Iterative CART Methods:
%   cart                - General interface for the CART methods
%   columnaction        - Column-action method with cyclic column sweeps
%
% Iterative SIRT Methods:
%   sirt                - General interface for all SIRT methods
%   cav                 - Component averaging (CAV) method
%   cimmino             - Cimmino's method
%   drop                - Diagonally Relaxed Orthogonal Projections (DROP) method
%   landweber           - Landweber's method
%   sart                - Simultaneous Algebraic Reconstruction Technique (SART) method
%
% Training Routines:
%   train_dpme          - Training method for the stopping rules DP and ME
%   train_relaxpar      - Training to find optimal relaxpar for ART/CART/SIRT methods
%
% Test Problems:
%   fancurvedtomo       - Creates 2D fan-beam curved-detector tomography test problem
%   fanlineartomo       - Creates 2D fan-beam linear-detector tomography test problem
%   paralleltomo        - Creates a 2D parallel-beam tomography test problem
%   phantomgallery      - A collection of 2D phantoms for use in test problems
%   purge_rows          - Remove zero or very sparse rows of A and corresp. entries in b
%   seismictomo         - Creates a 2D seismic travel-time tomography test problem
%   seismicwavetomo     - Seismic tomography problem without the ray assumption
%   show_tomo           - Illustrate the geometry of a tomographic test problem
%   sphericaltomo       - Creates a 2D spherical Radon tomography test problem
%
% Demonstration Scripts:
%   demo_art            - Demonstrates the use of, and the results from, the ART methods
%   demo_astra_2d       - Demonstrates interfacing to the ASTRA Toolbox
%   demo_cart           - Demonstrates the use of, and the results from, the CART method
%   demo_constraints    - Demonstrates the use of various constraints
%   demo_custom_all     - Demonstrates how to specify custom ART/SIRT methods
%   demo_custom_sirt    - Demonstrates the custom SIRT method
%   demo_matrixfree     - Demonstrates the matrix-free mode
%   demo_relaxpar       - Demonstrates the use of the relaxation parameter
%   demo_show_tomo      - Illustrates the use of the show_tomo function
%   demo_sirt           - Demonstrates the use of, and results from, SIRT methods
%   demo_stoprules      - Demonstrates how to use different stopping rules
%   demo_training       - Demonstrates the use of the training methods
%
% Auxiliary Routines:
%   afun_astra_2d_gpu   - Wrap ASTRA projectors into a function handle
%   afun_matrix         - Wrap a matrix into a function handle
%   calc_relaxpar       - Compute the relaxation parameter for all methods
%   check_inputs        - Check inputs and set default values
%   check_stoprules     - Check if stopping rule criteria are met
%   fbp                 - Filtered Back Projection using the system matrix
%   get_mfun_dfun       - Computes matrices M and D for SIRT methods