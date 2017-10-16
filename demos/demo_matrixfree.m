%DEMO_MATRIXFREE  Demonstrates the matrix-free mode
%
% This script illustrates the use of the default SIRT method sart applied
% to the same test problem using the matrix version, a "pseudo" matrix-free
% version in which the matrix exists but is wrapped using a function handle,
% and a fully matrix-free version in which the matrix is never explicitly
% stored.  The fully matrix-free mode is useful when storing the matrix
% explicitly requires too much memory.  Computations are instead carried
% out on the fly, i.e., in every application of forward/backprojection, so
% computing time is substantially longer.
%
% The script creates a parallel-beam CT test problem and solves it with
% the SART method in matrix, pseudo matrix-free and fully matrix-free mode.
% It further demonstrates the use of a function handle.  The exact phantom,
% the (identical) results of the methods and the (identical) sinograms and
% backprojected images are shown.
%
% See also: demo_astra_2d, demo_art, demo_cart, demo_constraints,
% demo_custom_all, demo_relaxpar, demo_sirt, demo_stoprules, demo_training.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_matrixfree:\n\n');

% Set the parameters for the test problem.
N = 50;           % The discretization points.
theta = 0:5:175;  % No. of used angles.
p = 75;           % No. of parallel rays.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.\n\n',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the matrix version of the test problem.
[A,b,x] = paralleltomo(N,theta,p);

% Show the exact solution
figure(1), clf
subplot(2,2,1)
imagesc(reshape(x,N,N)), colormap gray, axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 50;

fprintf(1,'Perform k = %2.0f iterations with the SART method.\n',k);

% Perform the SART iterations.
Xsart = sart(A,b,k);

% Show the SART solution.
subplot(2,2,2)
imagesc(reshape(Xsart,N,N)), colormap gray, axis image off
caxis(c)
title('SART')

% Pseudo matrix-free version by wrapping existing A into function handle
Afun_ps = @(XX,TT) afun_matrix(XX,TT,A);

fprintf(1,'Perform k = %2.0f iterations with the pseudo matrix-free SART method.\n',k);

% Perform the SART iterations.
Xsart_ps = sart(Afun_ps,b,k);

% Show the pseudo matrix-free SART solution.
subplot(2,2,3)
imagesc(reshape(Xsart_ps,N,N)), colormap gray, axis image off
caxis(c)
title('Pseudo matrix-free SART')

% Fully matrix-free version by matrix-free mode of paralleltomo.
% Create the matrx-version of the test problem.
Afun_mf = paralleltomo(N,theta,p,[],[],0);

fprintf(1,'Perform k = %2.0f iterations with the fully matrix-free SART method.\n',k);
fprintf(1,'  This is much slower - we trade speed for memory use!\n');

% Perform the SART iterations.
Xsart_mf = sart(Afun_mf,b,k);

% Show the pseudo matrix-free SART solution.
subplot(2,2,4)
imagesc(reshape(Xsart_mf,N,N)), colormap gray, axis image off
caxis(c)
title('Fully matrix-free SART')

%% Demonstrate usage of the matrix-free mode

% Forward projection, equivalent to A*x_ex
b_ps = Afun_ps(x, 'notransp');
b_mf = Afun_mf(x, 'notransp');

% Show the sinograms created by matrix, pseudo matrix-free and fully
% matrix-free forward projection.
figure(2), clf
subplot(1,3,1)
imagesc(reshape(b,p,length(theta)))
axis image off
title('Sinogram')
subplot(1,3,2)
imagesc(reshape(b_ps,p,length(theta)))
axis image off
title('Pseudo matrix-free')
subplot(1,3,3)
imagesc(reshape(b_mf,p,length(theta)))
axis image off
title('Fully matrix-free')

% Back-projection (unfiltered), equivalent to A'*b_ex.
z = A'*b;
z_ps = Afun_ps(b, 'transp');
z_mf = Afun_mf(b, 'transp');

% Show the images created by matrix, pseudo matrix-free and fully matrix-
% free back-projection.
figure(3), clf
subplot(1,3,1)
imagesc(reshape(z,N,N))
axis image off
title('Backprojection with A')
subplot(1,3,2)
imagesc(reshape(z_ps,N,N))
axis image off
title('Pseudo matrix-free')
subplot(1,3,3)
imagesc(reshape(z_mf,N,N))
axis image off
title('Fully matrix-free')

fprintf(1,'\n');
fprintf(1,'Compare sizes of matrix, pseudo and fully matrix-free operators.\n');

% Size of operator, equivalent to size(A).
fprintf(1,'Size of A matrix from the three different versions\n');
fprintf(1,'A is a matrix:\n');
disp(size(A))
fprintf(1,'A is a matrix wrapped into a function handle:\n');
disp(Afun_ps([], 'size'))
fprintf(1,'A is a true fundtion handle (no matrix):\n');
disp(Afun_mf([], 'size'))
