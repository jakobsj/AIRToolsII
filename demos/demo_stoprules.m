%DEMO_STOPRULES  Demonstrates how to use different stopping rules
%
% This script illustrates the use of two different stopping rules for
% Cimmino's method.
%
% The script creates a parallel-beam test problem, adds noise and solves
% the problem without two different stopping rules.  Then the error
% history is plotted together with the found number of iterations.
%
% See also: demo_art, demo_constraints, demo_custom_all, demo_matrixfree,
% demo_relaxpar, demo_sirt, demo_training, train_dpme.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_stoprules:\n\n');

% Problem specification.
N = 75;           % The image is N-by-N.
theta = 0:2:178;  % Projection angles.
p = 75;           % No. of parallel rays.
eta = 0.02;       % Relative noise level.
kmax = 1000;      % Number of iterations.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Generate test test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);
[m,n] = size(A);

% Add noise; we scale the noise vector e such that
% || e ||_2 / || b_ex ||_2 = eta.
rng(0);                         % Initialize random number generator.
e = randn(size(b_ex));          % Gaussian white noise.
e = eta*norm(b_ex)*e/norm(e);   % Scale the noise vector.
b = b_ex + e;                   % Add the noise to the pure data.

% Run Cimmino's method without a stopping rule, and compute the
% error history and its minimum.
fprintf(1,'\n\nRunning %2.0f iterations of Cimmino''s method',kmax);
X = cimmino(A,b,1:kmax);
E = zeros(kmax,1);
for k=1:kmax
    E(k) = norm(x_ex-X(:,k));
end
E = E/norm(x_ex);
[minE,K] = min(E);

% Run Cimmino's method with the discrepancy principle stopping rule.
% Since we know the exact norm of the noise, we use a tau close to 1.
fprintf(1,'\n\nRunning cimmino with DP stopping rule\n');
options.stoprule.type = 'DP';
options.stoprule.taudelta = 1.01*norm(e);
[~,info] = cimmino(A,b,1:kmax,[],options);
kDP = info.finaliter;
disp(['  DP found k = ',num2str(kDP)])

% Run Cimmino's method with the NCP stopping rule suited for 2D CT problem.
fprintf(1,'\nRunning cimmino with NCP stopping rule\n');
ntheta = length(theta);
p = m/ntheta;
options.stoprule.type = 'NCP';
options.stoprule.res_dims = [p,ntheta];
[~,info] = cimmino(A,b,1:kmax,[],options);
k2D = info.finaliter;
disp(['   NCP found k = ',num2str(k2D)])

figure(1), clf
plot(1:kmax,E,'-b',K,minE,'*b',kDP,E(kDP),'or',k2D,E(k2D),'dm',...
    'linewidth',1.5,'markersize',6)
legend('Relative error','Minimum','DP','NCP')
xlabel('Iterations')
ylabel('Relative error')