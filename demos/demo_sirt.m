%DEMO_SIRT  Demonstrates the use of, and results from, SIRT methods
%
% This script illustrates the use of the SIRT methods landweber,
% cimmino, cav, drop, and sart.
%
% The script creates a parallel-beam CT test problem and solves it with the
% SIRT methods. The exact solution and the results from the methods are shown.
%
% See also: demo_art, demo_constraints, demo_custom_all, demo_matrixfree,
% demo_relaxpar, demo_stoprules, demo_training.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_sirt:\n\n');

% Set the parameters for the test problem.
N = 50;           % The discretization points.
theta = 0:5:179;  % No. of used angles.
p = 75;           % No. of parallel rays.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the test problem.
[A,b,x] = paralleltomo(N,theta,p);

% Show the exact solution
figure(1), clf
subplot(2,3,1)
imagesc(reshape(x,N,N)), colormap gray, axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 50;

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Landweber''s method.\n',k);

% Perform the landweber iterations.
Xland = landweber(A,b,k);

% Show the landweber solution.
subplot(2,3,2)
imagesc(reshape(Xland,N,N)), colormap gray, axis image off
caxis(c)
title('Landweber')

fprintf(1,'Perform k = %2.0f iterations with Cimmino''s method.\n',k);

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k);

% Show the cimmino solution.
subplot(2,3,3)
imagesc(reshape(Xcimp,N,N)), colormap gray, axis image off
caxis(c)
title('Cimmino')

fprintf(1,'Perform k = %2.0f iterations with the CAV method.\n',k);

% Perform the CAV iterations.
Xcav = cav(A,b,k);

% Show the CAV solution.
subplot(2,3,4)
imagesc(reshape(Xcav,N,N)), colormap gray, axis image off
caxis(c)
title('CAV')

fprintf(1,'Perform k = %2.0f iterations with the DROP method.\n',k);

% Perform the DROP iterations.
Xdrop = drop(A,b,k);

% Show the DROP solution.
subplot(2,3,5)
imagesc(reshape(Xdrop,N,N)), colormap gray, axis image off
caxis(c)
title('DROP')

fprintf(1,'Perform k = %2.0f iterations with the SART method.\n',k);

% Perform the SART iterations.
Xsart = sart(A,b,k);

% Show the SART solution.
subplot(2,3,6)
imagesc(reshape(Xsart,N,N)), colormap gray, axis image off
caxis(c)
title('SART')