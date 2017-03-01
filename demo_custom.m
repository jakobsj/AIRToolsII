%DEMO_CUSTOM (script) Demonstrates how to specify custom ART and SIRT methods.
% See also: demo_art, demo_constraints, demo_sirt, demo_training.
% 
% The script creates a parallel-beam test problem, add noise and solves
% the problem with the the standard kaczmarz method and a custom ART method
% running through rows backwards. It also solves the problem using the
% standard landweber method as well as two custom variants by defining M
% and T matrices used in the general formulation of SIRT. The exact
% solution and the results from the methods are shown. 
%
% See also: demo_art, demo_constraints, demo_matrixfree, demo_sirt, 
% demo_training.

% Jakob Sauer Jorgensen, 2017-03-01, DTU Compute.

%% Set up test problem

close all
fprintf(1,'\nStarting demo_custom:\n\n');

% Set the parameters for the test problem.
N = 50;           % The discretization points.
theta = 0:5:179;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.05;       % Relative noise level.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the matrx-version of the test problem.
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

%% Standard kaczmarz

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a while ...');

% Perform the standard kaczmarz iterations using general interface
Xkacz = art('kaczmarz',A,b,k);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacz,N,N)), colormap gray,
axis image off
caxis(c);
title('Kaczmarz reconstruction')

%% Custom (backwards) kaczmarz

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations of custom (backwards) Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a while ...');

% Perform a custom Kaczmarz method running backwards through rows
row_order = size(A,1):-1:1;
Xkacz_back = art(row_order,A,b,k);

% Show the custom (backwards) kaczmarz solution.
figure
imagesc(reshape(Xkacz_back,N,N)), colormap gray,
axis image off
caxis(c);
title('Custom (backwards) Kaczmarz reconstruction')

%% Standard landweber

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with landweber method.',k);
fprintf(1,'\nThis takes a moment ...');

% Perform the standard landweber iteration using general interface
Xland = sirt('landweber',A,b,k);

% Show the landweber solution.
figure
imagesc(reshape(Xland,N,N)), colormap gray,
axis image off
caxis(c);
title('Landweber reconstruction')

%% Custom SIRT (landweber done manually to show interface)

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations of custom (also landweber) SIRT method.',k);
fprintf(1,'\nThis takes a moment ...');

% Set M and T function handles to do nothing as in landweber's method
Mfun = @(XX) XX;
Tfun = @(XX) XX;

% Perform the standard landweber iteration using general interface
Xland_c = sirt({Mfun,Tfun},A,b,k);

% Show the custom (backwards) kaczmarz solution.
figure
imagesc(reshape(Xland_c,N,N)), colormap gray,
axis image off
caxis(c);
title('Custom (also landweber) SIRT reconstruction')

%% Custom SIRT (landweber with arbitrary scaling)

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations of custom (scaled landweber) SIRT method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Set M and T function handles to do nothing as in landweber's method
scaling = 10;
Mfun = @(XX) scaling*XX;
Tfun = @(XX) (1/scaling)*XX;

% Perform the standard landweber iteration using general interface
Xland_s = sirt({Mfun,Tfun},A,b,k);

% Show the custom (backwards) kaczmarz solution.
figure
imagesc(reshape(Xland_s,N,N)), colormap gray,
axis image off
caxis(c);
title('Custom (scaled landweber) SIRT reconstruction')