%DEMO_CUSTOM_SIRT (script) Demonstrates the custom SIRT method
% 
% This script deomonstrates that the symmetric Kaczmarz algorithm is,
% theoretically, identical to a certain SIRT algorithm with a dense
% "weight" matrix M.
%
% The script creates a parallel-beam CT test problem and solves it with
% both the symmetric Kaczmarz method and the custom SIRT method with the
% specific matrix M.
%
% See also: demo_art, demo_constraints, demo_matrixfree, demo_sirt, 
% demo_training.

% PCH.

%% Set up test problem

clear, clc
rng(0);
figure(1), clf

% Set the parameters for the test problem.
N = 50;           % The image is N-by-N.
theta = 0:5:179;  % No. of used angles.
p = 75;           % No. of parallel rays.

% Create the test problem and purge zero rows from A.
[A,b,x] = paralleltomo(N,theta,p);
[A,b] = purge_rows(A,b);

% Show the exact solution.
subplot(2,2,1)
imagesc(reshape(x,N,N)), colormap gray, axis image off
title('Exact phantom','fontsize',14)

% No. of iterations.
kmax = 20;
omega = 1;

%% Symmetric Kaczmarz.

% Perform the symmetric Kaczmarz iterations.
options.relaxpar = omega;
Xsym = art('symkaczmarz',A,b,2:2:kmax,[],options);

% Show the Kaczmarz solution.
subplot(2,2,3)
imagesc(reshape(Xsym(:,end),N,N)), colormap gray, axis image off
caxis([0 1])
title('Symmetric Kaczmarz','fontsize',14)

%% Custom SIRT -> symmetric Kaczmarz.

% Compute the  "weight" matrix M that corresponds to symmetric Kaczmarz.
AAT = full(A*A');
L = tril(AAT,-1);
Delta = diag(diag(AAT));
M = (omega*Delta+L')\((2*omega-1)*Delta)/(omega*Delta+L);

% Set M and D function handles (D is the identity).
Mfun = @(XX)M*XX;
Dfun = @(XX) XX;

% Perform custom SIRT iterations using the general interface.  Note that
% one iteration here corresponds to a pair of iterations in sym. Kaczmarz.
Xsirt = sirt({Mfun,Dfun},A,b,1:1:kmax/2,[],options);

% Show the custom Landweber solution.
subplot(2,2,4)
imagesc(reshape(Xsirt(:,end),N,N)), colormap gray, axis image off
caxis([0 1]);
title('Custom SIRT','fontsize',14)

%% Compute and show the differences between the iteration vectors.

d = zeros(kmax/2,1);
for i=1:kmax/2;
    d(i) = norm(Xsym(:,i)-Xsirt(:,i))/norm(Xsym(:,i));
end
subplot(2,2,2)
plot(2:2:kmax,d,'o-','linewidth',1)
title('Differences')
xlabel('\itk')
set(gca,'fontsize',14)