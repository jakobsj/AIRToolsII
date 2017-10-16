%DEMO_ART Demonstrates the use of, and results from, the ART methods
%
% This script illustrates the use of the ART methods Kaczmarz, symmetric
% Kaczmarz, and randomized Kaczmarz.
%
% The script creates a parallel-beam CT test problem and solves it with the
% ART methods.  The exact image and the results from the methods are shown.
%
% See also: demo_cart, demo_constraints, demo_custom_all, demo_matrixfree,
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
fprintf(1,'Starting demo_art:\n\n');

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
imagesc(reshape(x,N,N)), colormap gray, axis image off
c = caxis;
title('Exact phantom')

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.\n',k);

% Perform the Kaczmarz iterations.
Xkacz = kaczmarz(A,b,k);

% Show the Kaczmarz solution.
subplot(2,2,2)
imagesc(reshape(Xkacz,N,N)), colormap gray, axis image off
caxis(c)
title('Kaczmarz')

fprintf(1,'Perform k = %2.0f iterations with the symmetric Kaczmarz method.\n',k);

% Perform the symmetric Kaczmarz iterations.
Xsymk = symkaczmarz(A,b,k);

% Show the symmetric kaczmarz solution.
subplot(2,2,3)
imagesc(reshape(Xsymk,N,N)), colormap gray, axis image off
caxis(c)
title('Symmetric Kaczmarz')

fprintf(1,'Perform k = %2.0f iterations with the randomized Kaczmarz method.\n',k);

% Perform the randomized Kaczmarz iterations.
Xrand = randkaczmarz(A,b,k);

% Show the randomized Kaczmarz solution.
subplot(2,2,4)
imagesc(reshape(Xrand,N,N)), colormap gray, axis image off
caxis(c)
title('Randomized Kaczmarz')
