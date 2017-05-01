%DEMO_CUSTOM_ALL (script) Demonstrates how to specify custom ART and SIRT methods
% 
% This script creates a parallel-beam CT test problem and solves the problem
% with the the standard Kaczmarz method and a custom ART method running
% through the rows backwards. It also solves the problem using the standard
% Landweber method and two custom variants by defining M and T matrices used
% in the general formulation of SIRT, the last of which makes the custom
% method identical to Cimmino (which we also show). The exact solution and
% the results from all the methods are shown.
%
% See also: demo_art, demo_constraints, demo_matrixfree, demo_sirt, 
% demo_training.

% Jakob Sauer Jorgensen, 2017-03-01, DTU Compute.

%% Set up test problem.

clear, clc
fprintf(1,'Starting demo_custom:\n\n');

% Set the parameters for the test problem.
N = 50;           % The image is N-times-N..
theta = 0:2:178;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.0;       % Relative noise level.

figure(1), clf

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the matrix version of the test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);

% Add noise to the rhs.
rng(0);
e = randn(size(b_ex));
b = b_ex + eta*norm(b_ex)*e/norm(e);

% No. of iterations, and color scale.
k_art = 10;
k_sirt = 50;
c = [0,1];

%% Standard Kaczmarz.

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.\n',k_art);

% Perform the standard Kaczmarz iterations using the general interface.
Xkacz = art('kaczmarz',A,b,k_art);

% Show the Kaczmarz solution.
subplot(2,4,1)
imagesc(reshape(Xkacz,N,N)), colormap gray, axis image off
caxis(c)
title('Kaczmarz')

%% Custom backwards Kaczmarz.

fprintf(1,'Perform k = %2.0f iterations of custom (backwards) Kaczmarz''s method.\n\n',k_art);

% Perform a custom Kaczmarz method running backwards through rows.
row_order = size(A,1):-1:1;
Xkacz_back = art(row_order,A,b,k_art);

% Show the custom (backwards) Kaczmarz solution.
subplot(2,4,5)
imagesc(reshape(Xkacz_back,N,N)), colormap gray, axis image off
caxis(c)
title('Custom (backwards) Kaczmarz')

%% Standard column-action method.

fprintf(1,'Perform k = %2.0f iterations with standard column-action method.\n',k_art);

% Perform the standard column-action iterations using the general interface
Xcolu = cart('columnaction',A,b,k_art);

% Show the column-action solution.
subplot(2,4,2)
imagesc(reshape(Xcolu,N,N)), colormap gray, axis image off
caxis(c)
title('Column action')

%% Custom (backwards) column action

fprintf(1,'Perform k = %2.0f iterations of custom (backwards) column-action method.\n\n',k_art);

% Perform a custom column-action method running backwards through rows
col_order = size(A,2):-1:1;
Xcolu_back = cart(col_order,A,b,k_art);

% Show the custom (backwards) column-action solution.
subplot(2,4,6)
imagesc(reshape(Xcolu_back,N,N)), colormap gray, axis image off
caxis(c)
title('Custom (backwards) column action')

%% Standard Landweber.

fprintf(1,'Perform k = %2.0f iterations with Landweber''s method.\n',k_sirt);

% Perform the standard Landweber iteration using the general interface.
Xland = sirt('landweber',A,b,k_sirt);

% Show the Landweber solution.
subplot(2,4,3)
imagesc(reshape(Xland,N,N)), colormap gray, axis image off
caxis(c)
title('Landweber')

%% Custom SIRT (Landweber done manually to illustrate the interface).

fprintf(1,'Perform k = %2.0f iterations of custom (= Landweber) SIRT method.\n',k_sirt);

% Set M and D function handles to identities.
Mfun = @(XX) XX;
Dfun = @(XX) XX;

% Perform the standard Landweber iteration using the general interface.
Xland_c = sirt({Mfun,Dfun},A,b,k_sirt);

% Show the custom Landweber solution.
subplot(2,4,7)
imagesc(reshape(Xland_c,N,N)), colormap gray, axis image off
caxis(c)
title('Unit weights -> Landweber')

%% Standard Cimmino.

fprintf(1,'Perform k = %2.0f iterations of Cimmino''s method.\n',k_sirt);

% Use the general interface.
Xland_s = sirt('cimmino',A,b,k_sirt);

% Show the Cimmino solution.
subplot(2,4,4)
imagesc(reshape(Xland_s,N,N)), colormap gray, axis image off
caxis(c)
title('Cimmino')

%% Custom SIRT (scaled Landweber -> Cimmino)

fprintf(1,'Perform k = %2.0f iterations of custom (scaled Landweber) SIRT method.\n\n',k_sirt);

% Set M and D function handles to achieve Cimmino's method.
M = 1./full(sum(A.^2,2));
Mfun = @(XX) M.*XX;
Dfun = @(XX) XX;

% Use the general interface.
Xland_s = sirt({Mfun,Dfun},A,b,k_sirt);

% Show the custom solution.
subplot(2,4,8)
imagesc(reshape(Xland_s,N,N)), colormap gray, axis image off
caxis(c)
title('Scaled Landweber -> Cimmino')