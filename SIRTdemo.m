%SIRTdemo (script) Demonstrates the use of, and the results from, the SIRT methods.
%
% This script illustrates the use of the SIRT methods landweber,
% cimmino, cav, drop, and sart.
%
% The script creates a parallel-beam test problem, add noise and solves
% the problem with the SIRT methods. The exact solution and the results
% from the methods are shown. 
%
% See also: ARTdemo, nonnegdemo, trainingdemo.

% Maria Saxild-Hansen and Per Chr. Hansen, May 23, 2010, DTU Compute.

close all
fprintf(1,'\nStarting SIRTdemo:\n\n');

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
randn('state',0);
e = randn(size(b_ex));
e = delta*e/norm(e);
b = b_ex + e;

% Show the exact solution
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 50;

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Landweber''s method.',k);
fprintf(1,'\nThis only takes a moment ...\n');

% Perform the landweber iterations.
Xland = landweber(A,b,k);

% Show the landweber solution.
figure
imagesc(reshape(Xland,N,N)), colormap gray,
axis image off
caxis(c);
title('Landweber reconstruction')

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with Cimmino''s method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the cimmino iterations.
Xcimp = cimmino(A,b,k);

% Show the cimmino solution.
figure
imagesc(reshape(Xcimp,N,N)), colormap gray,
axis image off
caxis(c);
title('Cimmino reconstruction')

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with the CAV method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the CAV iterations.
Xcav = cav(A,b,k);

% Show the CAV solution.
figure
imagesc(reshape(Xcav,N,N)), colormap gray,
axis image off
caxis(c);
title('CAV reconstruction')

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with the DROP method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the DROP iterations.
Xdrop = drop(A,b,k);

% Show the DROP solution.
figure
imagesc(reshape(Xdrop,N,N)), colormap gray,
axis image off
caxis(c);
title('DROP reconstruction')

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with the SART method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the SART iterations.
Xsart = sart(A,b,k);

% Show the SART solution.
figure
imagesc(reshape(Xsart,N,N)), colormap gray,
axis image off
caxis(c);
title('SART reconstruction')