%DEMO_CART (script) Demonstrates the use of, and the results from, the CART method.
%
% This script illustrates the use of the columnwise kaczmarz (CART) method 
% columnkaczmarz and for comparison a standard (rowwise) kaczmarz.
%
% The script creates a parallel-beam test problem, adds noise, and solves
% the problems with the ART methods.  The exact solution and the results
% from the methods are shown.
%
% See also: demo_constraints, demo_custom, demo_matrixfree, demo_sirt, 
% demo_training.

% Jakob Sauer Jorgensen, 2017-03-01, DTU Compute.

close all
fprintf(1,'\nStarting demo_cart:\n\n');

% Set the parameters for the test problem.
N = 50;           % The discretization points.
theta = 0:5:179;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.05;       % Relative noise level.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);

% Noise level.
delta = eta*norm(b_ex);

% Add noise to the rhs.
rng(0);
e = randn(size(b_ex));
e = delta*e/norm(e);
b = b_ex + e;

% Show the exact solution.
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 10;

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with the columnwise Kaczmarz method.',k);
fprintf(1,'\nThis takes a moment ...');

% Perform the columnwise Kaczmarz iterations.
Xcart = columnkaczmarz(A,b,k);

% Show the columnwise kaczmarz solution.
figure
imagesc(reshape(Xcart,N,N)), colormap gray,
axis image off
caxis(c);
title('Columnwise Kaczmarz reconstruction')

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the kaczmarz iterations.
Xkacz = kaczmarz(A,b,k);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacz,N,N)), colormap gray,
axis image off
caxis(c);
title('Kaczmarz reconstruction')

