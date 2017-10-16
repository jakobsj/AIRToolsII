%DEMO_TRAINING  Demonstrates the use of the training methods
%
% This script demonstrates the use of the functions train_relaxpar and
% train_dpme.  We train the SIRT method cimmino and the ART method kaczmarz.
% For the SIRT method the stopping rule ME is used, and for the ART method
% the stopping rule DP is used.
%
% See also: demo_art, demo_constraints, demo_custom_all, demo_matrixfree,
% demo_relaxpar, demo_sirt, demo_stoprules, train_dpme, train_relaxpar.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clc
fprintf(1,'Starting demo_training:\n\n');

% Set the parameters for the test problem:
N = 50;           % Discretization points.
theta = 0:5:179;  % No. of angles.
p = 75;           % No. of parallel rays.
eta = 0.03;       % Relative noise level.
kmaxSIRT = 100;   % Max number of SIRT iterations.
kmaxART = 50;     % Max number of ART iterations.

fprintf(1,'Creating a test problem with parallel tomography\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f and p = %2.0f.\n\n',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);

% Create the test problem and add noise.
[A,b_ex,x_ex] = paralleltomo(N,theta,p);
rng(0);
e = randn(size(b_ex));
e = eta*norm(b_ex)*e/norm(e);
b = b_ex + e;
delta = norm(e);

% Define the SIRT and ART methods.
SIRTmethod = @cimmino;
ARTmethod = @kaczmarz;

% Train the SIRT relaxation parameter.
fprintf(1,'Training the relaxation parameter for the SIRT method.\n');
relaxparSIRT = train_relaxpar(A,b,x_ex,SIRTmethod,kmaxSIRT);
fprintf(1,['  Found relaxpar = ',num2str(relaxparSIRT),'\n']);
optionsSIRT.relaxpar = relaxparSIRT;

% Train the ART relaxation parameter.
fprintf(1,'Training the relaxation parameter for the ART method.\n');
relaxparART = train_relaxpar(A,b,x_ex,ARTmethod,kmaxART);
fprintf(1,['  Found relaxpar = ',num2str(relaxparART),'\n\n']);
optionsART.relaxpar = relaxparART;

% Stopping rule for the SIRT and ART methods.
typeSIRT = 'ME';
typeART = 'DP';
% Number of samples.
s = 5;

% Train the paramter tau in the stopping rules.
fprintf(1,'Training the stopping parameter for the SIRT method.\n');
tauSIRT = ...
  train_dpme(A,b_ex,x_ex,SIRTmethod,typeSIRT,delta,s,kmaxSIRT,optionsSIRT);
fprintf(1,['  Found tau = ',num2str(tauSIRT),'\n']);
fprintf(1,'Training the stopping parameter for the ART method.\n');
tauART = ...
  train_dpme(A,b_ex,x_ex,ARTmethod,typeART,delta,s,kmaxART,optionsART);
fprintf(1,['  Found tau = ',num2str(tauART),'\n\n']);

% Set the stopping rules for the SIRT and ART methods.
optionsSIRT.stoprule.type = typeSIRT;
optionsART.stoprule.type = typeART;

% Set the taudelta parameter for the SIRT and ART methods.
optionsSIRT.stoprule.taudelta = tauSIRT*delta;
optionsART.stoprule.taudelta = tauART*delta;
% optionsSIRT.stoprule.taudelta = 2*delta;
% optionsART.stoprule.taudelta = 2*delta;

% Iterating with the SIRT method.
fprintf(1,'Use the SIRT method with the found parameters.\n');
[XSIRT,infoSIRT] = SIRTmethod(A,b,kmaxSIRT,[],optionsSIRT);

% Iterating with the ART method.
fprintf(1,'Use the ART method with the found parameters.\n');
[XART,infoART] = ARTmethod(A,b,kmaxART,[],optionsART);

% Show the results.
figure(1), clf
subplot(2,2,1)
imagesc(reshape(x_ex,N,N)), colormap gray, axis image off
c = caxis;
title('Exact phantom')

subplot(2,2,3)
imagesc(reshape(XSIRT,N,N)), colormap gray, axis image off
caxis(c);
title(['SIRT: k = ',num2str(infoSIRT.finaliter)])

subplot(2,2,4)
imagesc(reshape(XART,N,N)), colormap gray, axis image off
caxis(c);
title(['ART: k = ',num2str(infoART.finaliter)])