%DEMO_CONSTRAINTS  Demonstrates the use of various constraints
%
% This script illustrates the use of Cimmino's method with the lbound
% and ubound options.
%
% The script creates a parallel-beam test problem, adds noise and solves
% the problem without and with the constraints.  The exact solution and
% the results from the methods are shown. The script shows that we obtain
% better reconstructions when we are able to impose bounds.
%
% See also: demo_art, demo_cart, demo_custom_all, demo_matrixfree,
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
fprintf(1,'Starting demo_constraints:\n\n');

% Set the parameters for the test problem.
N = 50;           % The image is N-by-N.
theta = 0:2:178;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.02;       % Relative noise level.
k = 5000;         % Number of iterations.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);

% Add noise to the data.
rng(0);
e = randn(size(b_ex));
b = b_ex + eta*norm(b_ex)*e/norm(e);

% Show the exact solution
figure(1), clf
subplot(2,3,1)
imagesc(reshape(x_ex,N,N)), colormap gray, axis image off
c = caxis;
title('Exact phantom')

% Perform Cimmino without constraints.
fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Cimmino''s method.\n',k);
Xcim = cimmino(A,b,k);

% Show the Cimmino solution.
subplot(2,3,2)
imagesc(reshape(Xcim,N,N)), colormap gray, axis image off
caxis(c)
title('Cimmino')

% Perform Cimmino iterations with nonnegativity constraints.
fprintf(1,'Perform k = %2.0f iterations with nonnegativity constraints.\n',k);
options.lbound = 0;
Xcimn = cimmino(A,b,k,[],options);

% Show the solution.
subplot(2,3,4)
imagesc(reshape(Xcimn,N,N)), colormap gray, axis image off
caxis(c)
title('With nonneg.')

% Perform Cimmino iterations with box constraints.
fprintf(1,'Perform k = %2.0f iterations with box constraints.\n',k);
options.lbound = 0;
options.ubound = 1;
Xcimb = cimmino(A,b,k,[],options);

% Show the solution.
subplot(2,3,5)
imagesc(reshape(Xcimb,N,N)), colormap gray,
axis image off
caxis(c)
title('With box')

% Determine the indices to those regions in the image that are known to
% have exact pixel values equal to 0.3, and enforce very tight bounds
% (simulation equality constraints) on these elements.
I = find(abs(x_ex-0.3)<1e-10);
L = zeros(N^2,1);
L(I) = 0.299;
U = ones(N^2,1);
U(I) = 0.301;
options.lbound = L;
options.ubound = U;

% Perform the cimmino iterations.
fprintf(1,'Repeat Cimmino with tight lower and upper bounds in certain pixels.\n\n');
Xcime = cimmino(A,b,k,[],options);

% Show the solution.
subplot(2,3,6)
imagesc(reshape(Xcime,N,N)), colormap gray, axis image off
caxis(c)
title('With box & tight')

% Show that enforcing constraints improves the reconstruction, by computing
% the errors in those pixels where the tight constraints are not active
z_ex = x_ex; z_ex(I) = [];
Zcim = Xcim; Zcim(I) = [];
Zcimn = Xcimn; Zcimn(I) = [];
Zcimb = Xcimb; Zcimb(I) = [];
Zcime = Xcime; Zcime(I) = [];
fprintf(1,'Errors in the reconstructions.\n')
fprintf(1,'First line over all pixels.\n')
fprintf(1,'Second line over those pixels where eq. constr. is not active\n')
fprintf(1,'No constraints, nonnegativity, box, and box + equality:\n');
disp([norm(x_ex-Xcim),norm(x_ex-Xcimn),norm(x_ex-Xcimb),norm(x_ex-Xcime)])
disp([norm(z_ex-Zcim),norm(z_ex-Zcimn),norm(z_ex-Zcimb),norm(z_ex-Zcime)])