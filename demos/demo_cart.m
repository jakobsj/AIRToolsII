%DEMO_CART  Demonstrates the use of, and results from, the CART method
%
% This script illustrates the use of the column-action reconstruction
% method and compares it with the standard (row-action) Kaczmarz method.
%
% The script creates a parallel-beam CT test problem and solves it with the
% CART and ART methods (for comparison).  The exact image and the results
% from the methods are shown.
%
% See also: demo_art, demo_constraints, demo_custom_all, demo_matrixfree,
% demo_relaxpar, demo_sirt, demo_stoprules, demo_training.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_cart:\n\n');

% Set the parameters for the test problem.
N = 50;           % The image is N-times-N..
theta = 0:2:178;  % No. of used angles.
p = 75;           % No. of parallel rays.
k = 10;           % Number of iterations.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the test problem.
[A,b,x] = paralleltomo(N,theta,p);

% Show the exact solution.
figure(1), clf
subplot(2,2,1)
imagesc(reshape(x,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with the column-action method.\n',k);

% Perform the column-action iterations.
Xcart = columnaction(A,b,k);

% Show the solution.
subplot(2,2,3)
imagesc(reshape(Xcart,N,N)), colormap gray, axis image off
caxis(c)
title('Column-action method')

fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.\n',k);

% Perform the Kaczmarz iterations.
Xkacz = kaczmarz(A,b,k);

% Show the solution.
subplot(2,2,4)
imagesc(reshape(Xkacz,N,N)), colormap gray, axis image off
caxis(c)
title('Kaczmarz')