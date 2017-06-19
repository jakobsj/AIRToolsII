%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIR Tools II
% Version 2.0 of AIR Tools
% June 2017
% 
% Per Christian Hansen and Jakob Sauer Jørgensen
% Dept. Applied Mathematics and Computer Science (DTU Compute)
% Technical University of Denmark
%
% Contact: pcha@dtu.dk
% 
% This file is part of the AIR Tools package and is distributed under
% the 3-Clause BSD Licence.  A separate license file should be provided
% as part of the package. 
% 
% Copyright 2017 Per Christian Hansen & Jakob Sauer Jørgensen, DTU Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Info:
%   CONTENTS.m          - This file with overview of available files.
%   CHANGES.txt         - List of changes made in code since last version.
%   LICENSE.txt         - License file.
%
% Iterative ART Methods:
%   art                 - General interface for all Kaczmarz/ART methods.
%   kaczmarz            - Kaczmarz's method with cyclic row sweep
%                         (often referred to as ART).
%   randkaczmarz        - Kaczmarz's method with random row selection.
%   symkaczmarz         - Kaczmarz's method with ``symmetric'' (top to
%                         bottom to top) row sweep.
%
% Iterative CART Methods:
%   cart                - General interface for the CART methods.
%   columnaction        - Column-action method with cyclic column sweeps.
%
% Iterative SIRT Methods:
%   sirt                - General interface for all SIRT methods.
%   cav                 - Component averaging (CAV) method.
%   cimmino             - Cimmino's method.
%   drop                - Diagonally relaxed orthogonal projections (DROP)
%                         method.
%   landweber           - Landweber's method.
%   sart                - Simultaneous algebraic reconstruction technique
%                         (SART) method.
%
% Training Routines:
%   train_dpme          - Use training to compute a good multiplicative
%                         factor in the DP and ME stopping rules.
%   train_relaxpar      - Use training to compute a good relaxation
%                         parameter for the ART/CART/SIRT methods.
%
% Test Problems:
%   fancurvedtomo       - 2D fan-beam CT with a curved detector (equal
%                         angles between rays).
%   fanlineartomo       - 2D fan-beam CT with a linear detector.
%   paralleltomo        - 2D parallel-beam CT.
%   phantomgallery      - A collection of different 2D phantoms.
%   purge_rows          - Remove empty or very sparse rows of A and the
%                         corresponding elements of b.
%   seismictomo         - 2D seismic travel-time tomography test problem.
%   seismicwavetomo     - Similar to seismictomo but without a ray assumption.
%   show_tomo           - Illustrate the geometry of a tomographic test
%                         problem using the rows of the matrix A.
%   sphericaltomo       - 2D spherical Radon transform tomography problem.
%
% Demonstration Scripts:
%   demo_art            - Demonstrates the use of, and the results from,
%                         the ART methods.
%   demo_astra_2d       - Demonstrates interfacing to the ASTRA Toolbox.
%   demo_cart           - Demonstrates the use of, and the results from,
%                         the CART method.
%   demo_constraints    - Demonstrates the use of nonnegativity and other
%                         constraints.
%   demo_custom_all     - Demonstrates how to specify custom ART and SIRT
%                         methods.
%   demo_custom_sirt    - Demonstrates how the custom SIRT interface can be
%                         used to implement the symmetric Kaczmarz method.
%   demo_matrixfree     - Demonstrates the matrix-free mode for test problems
%                         and reconstruction methods.
%   demo_relaxpar       - Demonstrates how to use various relaxation 
%                         parameter selection strategies.
%   demo_show_tomo      - Demonstrates tomography test problems.
%   demo_sirt           - Demonstrates the use of, and the results from,
%                         the SIRT methods.
%   demo_stoprules      - Demonstrates how to use different stopping rules.
%   demo_training       - Demonstrates the use of the training methods.
%
% Auxiliary Routines:
%   afun_astra_2d_gpu   - Wrap ASTRA projectors into a function handle.
%   afun_matrix         - Wrap a matrix into a function handle.
%   calc_relaxpar       - Compute the relaxation parameter for all methods.
%   check_inputs        - Check inputs and set default values.
%   check_stoprules     - Check if stopping rule criteria are met.
%   fbp                 - Filtered back projection using the system matrix.
%   get_mfun_dfun       - Computes matrices M and D for SIRT methods.