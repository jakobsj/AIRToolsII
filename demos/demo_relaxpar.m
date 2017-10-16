%DEMO_RELAXPAR Demonstrates the use of the relaxation parameter
%
% This script illustrates the use of different relaxation parameter
% strategies in Cimmino's method, including the one found by training.
%
% The script creates a parallel-beam test problem, adds noise and solves
% the problem with different choices of relaxation parameters.  Then the
% error histories are plotted.
%
% The script also illustrates that for noisy data, the line-search strategy
% leads to some oscillations in the error history.  These oscillations
% increase with the noise level.
%
% See also: train_relaxpar.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_relaxpar:\n\n');

% Problem specification.
N = 50;           % The image is N-by-N.
theta = 0:3:179;  % Projection angles.
p = 75;           % No. of parallel rays.
eta = 0.004;      % Relative noise level.
kmax = 500;       % Number of iterations.

fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.\n\n',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Generate test test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);

% Add noise; we scale the noise vector e such that
% || e ||_2 / || b_ex ||_2 = eta.
rng(0);                         % Initialize random number generator.
e = randn(size(b_ex));          % Gaussian white noise.
e = eta*norm(b_ex)*e/norm(e);   % Scale the noise vector.
b = b_ex + e;                   % Add the noise to the pure data.

% Get the largest allowed relaxation parameter and the one found by
% training.
fprintf(1,'Compute the largest allowed relaxation parameter\n');
[~,info] = cimmino(A,b,1);
relaxpar_max = 2/info.rho;
fprintf(1,'   relaxpar_max = %2.0f\n\n',relaxpar_max);
fprintf(1,'Compute the relaxation parameter via training\n');
relaxpar_train = train_relaxpar(A,b,x_ex,@cimmino,kmax);
fprintf(1,'   relaxpar_train = %2.0f\n\n',relaxpar_train);

% Try different strategies.
relaxpar = 0.5*relaxpar_max;
fprintf(1,'Running cimmino with fixed relaxpar = %2.0f\n',relaxpar);
options.relaxpar = relaxpar;
X = cimmino(A,b,1:kmax,[],options);
fprintf(1,'Running cimmino with trained relaxpar = %2.0f\n',relaxpar_train);
options.relaxpar = relaxpar_train;
Xtrain = cimmino(A,b,1:kmax,[],options);
fprintf(1,'Running cimmino with line strategy\n');
options.relaxpar = 'line';
Xline = cimmino(A,b,1:kmax,[],options);
options.relaxpar = 'psi2mod';
fprintf(1,['Running cimmino with ',options.relaxpar,' strategy\n']);
Xpsi1 = cimmino(A,b,1:kmax,[],options);
E = zeros(kmax,4);
for k=1:kmax
    E(k,1) = norm(x_ex-X(:,k));
    E(k,2) = norm(x_ex-Xtrain(:,k));
    E(k,3) = norm(x_ex-Xline(:,k));
    E(k,4) = norm(x_ex-Xpsi1(:,k));
end
E = E/norm(x_ex);

figure(1), clf
plot(1:kmax,E,'linewidth',1.5)
xlabel('Iterations')
ylabel('Relative error')
legend(['\omega = ',num2str(relaxpar)],...
       ['\omega_{train} = ',num2str(relaxpar_train)],...
       'line-search strategy',[options.relaxpar,' strategy'])   
% Repeat the experiment with 10 times larger noise.
b = b_ex + 10*e;
kmax = 50;

% Get the largest allowed relaxation parameter and the one found by
% training.
fprintf(1,'\n\nRepeat experiment with 10 times lareger noise\n');
fprintf(1,'Compute the largest allowed relaxation parameter\n');
[~,info] = cimmino(A,b,1);
relaxpar_max = 2/info.rho;
fprintf(1,'   relaxpar_max = %2.0f\n\n',relaxpar_max);
fprintf(1,'Compute the relaxation parameter via training\n');
relaxpar_train = train_relaxpar(A,b,x_ex,@cimmino,kmax);
fprintf(1,'   relaxpar_train = %2.0f\n\n',relaxpar_train);

% Try different strategies.
relaxpar = 0.5*relaxpar_max;
fprintf(1,'Running cimmino with fixed relaxpar = %2.0f\n',relaxpar);
options.relaxpar = relaxpar;
X = cimmino(A,b,1:kmax,[],options);
fprintf(1,'Running cimmino with trained relaxpar = %2.0f\n',relaxpar_train);
options.relaxpar = relaxpar_train;
Xtrain = cimmino(A,b,1:kmax,[],options);
fprintf(1,'Running cimmino with line strategy\n');
options.relaxpar = 'line';
Xline = cimmino(A,b,1:kmax,[],options);
options.relaxpar = 'psi1';
fprintf(1,['Running cimmino with ',options.relaxpar,' strategy\n']);
options.relaxpar = 'psi1';
Xpsi1 = cimmino(A,b,1:kmax,[],options);
E = zeros(kmax,4);
for k=1:kmax
    E(k,1) = norm(x_ex-X(:,k));
    E(k,2) = norm(x_ex-Xtrain(:,k));
    E(k,3) = norm(x_ex-Xline(:,k));
    E(k,4) = norm(x_ex-Xpsi1(:,k));
end
E = E/norm(x_ex);
[minE,K] = min(E);

fprintf(1,'\nNote the large oscillations in the error for ''line'' here\n');

figure(2), clf
plot(1:kmax,E,'linewidth',1.5)
xlabel('Iterations')
ylabel('Relative error')
legend(['\omega = ',num2str(relaxpar)],...
       ['\omega_{train} = ',num2str(relaxpar_train)],...
       'line-search strategy',[options.relaxpar,' strategy'])