%trainingdemo (script) Demonstrates the use of the training methods.
%
% This script demonstrates the use of the training functions
% trainLambdaSIRT, trainLambdaART, and trainDPME.  We train the SIRT
% method cimmino and the ART method kaczmarz.  For the SIRT method the
% stopping rule ME is used, and for the ART method the stopping rule DP
% is used.  Note that training the relaxation parameter for ART, training
% the stopping parameter for ART, and using ART takes several minutes.
%
% See also: ARTdemo, nonnegdemo, SIRTdemo.

% Maria Saxild-Hansen and Per Chr. Hansen, May 23, 2010, DTU Compute.

close all
fprintf(1,'\nStarting trainingdemo:\n\n');

% Set the parameters for the test problem:
N = 50;           % Discretization points.
theta = 0:5:179;  % No. of angles.
p = 75;           % No. of parallel rays.
eta = 0.05;       % Relative noise level.

fprintf(1,'Creating a test problem with parallel tomography\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f and p = %2.0f.\n',...
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

% Define the SIRT and ART methods.
SIRTmethod = @cimmino;
ARTmethod = @kaczmarz;

fprintf(1,'\nTraining the relaxation parameter for the SIRT method.');
fprintf(1,'\nThis only takes some seconds\n');

% Train the SIRT relaxation parameter.
lambdaSIRT = trainLambdaSIRT(A,b,x_ex,SIRTmethod);

% Set the relaxation parameter for the SIRT options.
optionsSIRT.lambda = lambdaSIRT;

fprintf(1,'\nTraining the relaxation parameter for the ART method.');
fprintf(1,'\nThis takes several minutes\n');

% Train the ART relaxation parameter.
lambdaART = trainLambdaART(A,b,x_ex,ARTmethod);

% Set the relaxation parameter for the SIRT options.
optionsART.lambda = lambdaART;

% Stopping rule for the SIRT and ART methods.
typeSIRT = 'ME';
typeART = 'DP';
% No of samples.
s = 5;

fprintf(1,'\nTraining the stopping parameter for the SIRT method.');
fprintf(1,'\nThis takes a few seconds\n');
tauSIRT = trainDPME(A,b_ex,x_ex,SIRTmethod,typeSIRT,delta,s,optionsSIRT);
fprintf(1,'\nTraining the stopping parameter for the ART method.');
fprintf(1,'\nThis takes several minutes\n');
tauART = trainDPME(A,b_ex,x_ex,ARTmethod,typeART,delta,s,optionsART);

% Set the stopping rules for the SIRT and ART methods.
optionsSIRT.stoprule.type = typeSIRT;
optionsART.stoprule.type = typeART;

% Set the taudelta parameter for the SIRT and ART methods.
optionsSIRT.stoprule.taudelta = tauSIRT*delta;
optionsART.stoprule.taudelta = tauART*delta;

fprintf(1,'\nUse the SIRT method.\n');
% Iterating with the SIRT method.
[XSIRT,infoSIRT] = SIRTmethod(A,b,1000,[],optionsSIRT);

fprintf(1,'\nUse the ART method.');
fprintf(1,'\nThis takes several minuts\n');
% Iterating with the ART method.
[XART,infoART] = ARTmethod(A,b,100,[],optionsART);

% Show the results.
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

figure
imagesc(reshape(XSIRT,N,N)), colormap gray,
axis image off
caxis(c);
title(['SIRT: k = ',num2str(infoSIRT(2))])

figure
imagesc(reshape(XART,N,N)), colormap gray,
axis image off
caxis(c);
title(['ART: k = ',num2str(infoART(2))])