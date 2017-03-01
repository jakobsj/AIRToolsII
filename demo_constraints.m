%DEMO_CONSTRAINTS (script) Demonstrates the use of nonnegativity and more
%general upper and lower bound constraints.
%
% This script illustrates the use of the cimmino and kaczmarz methods
% with the lbound and ubound options.
%
% The script creates a parallel-beam test problem, adds noise and solves
% the problem without and with the constraints.  The exact solution and
% the results from the methods are shown. 
%
% See also: demo_art, demo_custom, demo_matrixfree, demo_sirt, 
% demo_training.

% Maria Saxild-Hansen and Per Chr. Hansen, Oct. 21, 2010, DTU Compute.

close all
clear options
fprintf(1,'\nStarting demo_constraints:\n\n');

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

% Show the exact solution
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

% No. of iterations; set nonnegativity off. This is also the default
% behavior obtained if not specifying an lbound field of options.
k = 50;
options.lbound = -inf;

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Cimmino''s method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k,[],options);

% Show the cimmino solution.
figure
imagesc(reshape(Xcimp,N,N)), colormap gray,
axis image off
caxis(c);
title('Cimmino reconstruction')

% Set nonnegativity on.
options.lbound = 0;

fprintf(1,'\n');
fprintf(1,'Repeat with nonnegativity constraints.\n');
fprintf(1,'This takes a moment ...\n');

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k,[],options);

% Show the cimmino solution.
figure
imagesc(reshape(Xcimp,N,N)), colormap gray,
axis image off
caxis(c);
title('Cimmino reconstruction w. nonnegativity')

% Set nonnegativity off.
options.lbound = -inf;

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a while ...\n');

% Perform the kaczmarz iterations.
Xkacp = kaczmarz(A,b,k,[],options);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacp,N,N)), colormap gray,
axis image off
caxis(c);
title('Kaczmarz reconstruction')

% Set nonnegativity on.
options.lbound = 0;

fprintf(1,'\n');
fprintf(1,'Repeat with nonnegativity constraints.\n');
fprintf(1,'This takes a while ...\n');

% Perform the kaczmarz iterations.
Xkacp = kaczmarz(A,b,k,[],options);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacp,N,N)), colormap gray,
axis image off
caxis(c);
title('Kaczmarz reconstruction w. nonnegativity')

% Set scalar lower and upper bounds (too tight to demonstrate effect).
options.lbound = 0.1;
options.ubound = 0.7;

fprintf(1,'\n');
fprintf(1,'Repeat Cimmino with too tight lower and upper scalar bounds.\n');
fprintf(1,'This takes a moment ...\n');

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k,[],options);

% Show the cimmino solution.
figure
imagesc(reshape(Xcimp,N,N)), colormap gray,
axis image off
caxis(c);
title('Cimmino reconstruction w. bounds 0.1 and 0.7')

% Set vector lower and upper bounds: Different bounds in left and right
% half of image (all too tight to demonstrate effect). If non-scalar,
% vector input must have the same length as image to be reconstructed to
% enforce elementwise bounds.
ov = ones(N^2/2,1);
options.lbound = [0.1*ov; 0.15*ov];
options.ubound = [0.6*ov; 0.5*ov];

fprintf(1,'\n');
fprintf(1,'Repeat Cimmino with too tight lower and upper scalar bounds.\n');
fprintf(1,'This takes a moment ...\n');

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k,[],options);

% Show the cimmino solution.
figure
imagesc(reshape(Xcimp,N,N)), colormap gray,
axis image off
caxis(c);
title('Cimmino reconstruction w. vector bounds.')
