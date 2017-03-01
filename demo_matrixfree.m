%DEMO_MATRIXFREE (script) Demonstrates the matrix-free mode for test
%problems and reconstruction methods.
%
% As an example, this script illustrates the use of the default SIRT 
% method sart applied to the same test problem using the matrix version, a
% "pseudo" matrix-free version in which the matrix exists but is wrapped
% using a function handle, and a fully matrix-free version in which the
% matrix is never explicitly formed. The fully matrix-free mode can be
% useful when storing the full operator explicitly requires too much
% memory. Computations are instead carried out on the fly, i.e., in every
% application of forward/backprojection, so computation time will be
% substantially longer.
%
% The script creates a parallel-beam test problem, add noise and solves
% the problem with the sart method in matrix, pseudo matrix-free and fully
% matrix-free mode. It further demonstrates usage of the matrix-free
% operator represented using a function handle. The exact phantom, the
% (identical) results of the methods and the (identical) sinograms and
% backprojected images are shown.
%
% See also: demo_art, demo_constraints, demo_custom, demo_sirt, 
% demo_training.

% Jakob Sauer Jorgensen, 2017-03-01, DTU Compute.

close all
fprintf(1,'\nStarting demo_matrixfree:\n\n');

% Set the parameters for the test problem.
N = 50;           % The discretization points.
theta = 0:5:179;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.05;       % Relative noise level.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.\n',...
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

% Show the exact solution
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 50;

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

%% Pseudo matrix-free version by wrapping existing A into function handle
Afun_ps = @(XX,TT) afun_matrix(XX,TT,A);

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with the pseudo matrix-free SART method.',k);
fprintf(1,'\nThis takes a moment ...\n');

% Perform the SART iterations.
Xsart_pseudo = sart(Afun_ps,b,k);

% Show the pseudo matrix-free SART solution.
figure
imagesc(reshape(Xsart_pseudo,N,N)), colormap gray,
axis image off
caxis(c);
title('Pseudo matrix-free SART reconstruction')

%% Fully matrix-free version by matrix-free mode of paralleltomo.
% Create the matrx-version of the test problem.
is_matrix = 0;
Afun_mf = paralleltomo(N,theta,p,[],[],0);

fprintf(1,'\n');
fprintf(1,'Perform k = %2.0f iterations with the fully matrix-free SART method.',k);
fprintf(1,'\nThis takes a while ...\n');

% Perform the SART iterations.
Xsart_mf = sart(Afun_mf,b,k);

% Show the pseudo matrix-free SART solution.
figure
imagesc(reshape(Xsart_mf,N,N)), colormap gray,
axis image off
caxis(c);
title('Fully matrix-free SART reconstruction')

%% Demonstrate usage of the matrix-free mode

% Forward projection, equivalent to A*x_ex
b_ex_ps = Afun_ps(x_ex, 'notransp');
b_ex_mf = Afun_mf(x_ex, 'notransp');

% Show the sinograms created by matrix, pseudo matrix-free and fully
% matrix-free forward projection.
figure
subplot(1,3,1)
imagesc(reshape(b_ex,p,length(theta)))
axis image off
title('Sinogram matrix')
subplot(1,3,2)
imagesc(reshape(b_ex_ps,p,length(theta)))
axis image off
title('Pseudo matrix-free')
subplot(1,3,3)
imagesc(reshape(b_ex_mf,p,length(theta)))
axis image off
title('Fully matrix-free')

% Back-projection (unfiltered), equivalent to A^T*b_ex
z_ex = A'*b_ex;
z_ex_ps = Afun_ps(b_ex, 'transp');
z_ex_mf = Afun_mf(b_ex, 'transp');

% Show the images created by matrix, pseudo matrix-free and fully matrix-
% free back-projection.
figure
subplot(1,3,1)
imagesc(reshape(z_ex,N,N))
axis image off
title('Backprojection matrix')
subplot(1,3,2)
imagesc(reshape(z_ex_ps,N,N))
axis image off
title('Pseudo matrix-free')
subplot(1,3,3)
imagesc(reshape(z_ex_mf,N,N))
axis image off
title('Fully matrix-free')

fprintf(1,'\n');
fprintf(1,'Compare sizes of matrix, pseudo and fully matrix-free operators.\n');

% Size of operator, equivalent to size(A)
sizeA = size(A)
sizeA_ps = Afun_ps([], 'size')
sizeA_mf = Afun_mf([], 'size')
