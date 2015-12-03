% AIR Tools.
% Version 1.3  6-July-15
%
% Iterative ART Methods.
%   cart            - Columnwise version of Kaczmarz's method.
%   kaczmarz        - Kaczmarz's method (often referred to as ART).
%   randkaczmarz    - Randomized Kaczmarz method.
%   symkaczmarz     - Symmetric Kaczmarz method.
%
% Iterative SIRT Methods.
%   cav             - Component Averaging (CAV) method.
%   cimmino         - Cimmino's method.
%   drop            - Diagonally Relaxed Orthogonal Projections (DROP) method.
%   landweber       - The classical Landweber method.
%   sart            - Simultaneous Algebraic Reconstruction Technique (SART) method.
%
% Training Routines.
%   trainDPME       - Training method for the stopping rules DP and ME.
%   trainLambdaART  - Training to determine optimal lambda for ART methods.
%   trainLambdaSIRT - Training to determine optimal lambda for SIRT method.
%
% Test Problems.
%   fanbeamtomo     - 2D tomography test problem using fan beams.
%   paralleltomo    - 2D tomography test problem using parallel beams.
%   phantomgallery  - A collection of different 2D phantoms.
%   seismictomo     - 2D seismic travel-time tomography test problem.
%   seismicwavetomo - Similar to seismictomo but without a ray assumption.
%
% Demonstration Scripts.
%   ARTdemo         - Demonstrates the use of, and the results from, the ART methods.
%   nonnegdemo      - Demonstrates the use of nonnegativity constraints.
%   SIRTdemo        - Demonstrates the use of, and the results from, the SIRT methods.
%   trainingdemo    - Demonstrates the use of the training methods.
%
% Auxiliary Routines.
%   calczeta        - Calculates a specific root of a certain polynomial.
%   fbp             - Similar to iradon, but conforms with AIR Tools variables.
%   rzr             - Remove zero rows of A and the corresponding elements of b.