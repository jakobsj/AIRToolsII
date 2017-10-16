%DEMO_ASTRA_2D  Demonstrates interfacing to the ASTRA Toolbox
%
% This script demonstrates how to interface from AIR Tools to other
% software, in the present case to the ASTRA Tomography Toolbox, see 
% http://www.astra-toolbox.com/. Instead of using the matrices or matrix-
% free versions of the test problems in AIR Tools, ASTRA can be used for
% efficient GPU-accelerated application of the matrix and its transpose
% from all AIR Tools reconstruction functions.
%
% This script compares a fanlineartomo case reconstructed by the SART
% algorithm using the matrix and matrix-free modes of AIR Tools,
% interfacing from AIR Tools to ASTRA for application of the matrix and
% its transpose, as well as the built-in "SIRT" algorithm of ASTRA which is
% equivalent to AIR Tools' SART.
%
% See also: afun_astra_2d, demo_sirt.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_astra_2d:\n\n');

% Parameters for fanlineartomo test problem
N = 128;
dtheta = 1;
theta = 0:dtheta:(360-dtheta);
p = 1.5*N;
R = 2;
dw = 2;
sd = 3;

% Test image
X = phantomgallery('shepplogan',N);
x = X(:);

% AIRtools matrix and matrix-free
fprintf(1,'Creating a fanlineartomo tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);
A = fanlineartomo(N,theta,p,R,dw,sd);
Afun = fanlineartomo(N,theta,p,R,dw,sd,[],0);

% Data
b = A*x;
B = reshape(b,p,length(theta));

% Set up function calling ASTRA forward and backprojectors.
proj_geom = astra_create_proj_geom('fanflat', dw*N/p, p, ...
    (pi/180)*theta,R*N,(sd-R)*N);
vol_geom = astra_create_vol_geom(N,N);
astra_fun_gpu = @(XX,TT) afun_astra_2d_gpu(XX,TT,proj_geom,vol_geom);

% AIR Tools options to print progress, choice of relaxation parameter and
% number of SART iterations to run
opts.verbose = 1;
opts.relaxpar = 1;

num_iters = 10;

% AIR Tools SART with matrix
fprintf(1,'\n');
fprintf(1,...
    'Perform k = %2.0f iterations of AIR Tools sart with matrix.\n',...
    num_iters);
[Xsart_air,info_air] = sart(A,b,1:num_iters,[],opts);

% AIR Tools SART matrix-free
fprintf(1,'\n');
fprintf(1,['Perform k = %2.0f iterations of AIR Tools sart ',...
    'matrix-free (takes a while...).\n'],...
    num_iters);
[Xsart_air_mf,info_air_mf] = sart(Afun,b,1:num_iters,[],opts);

% AIR Tools SART with ASTRA GPU matrix-free operator
fprintf(1,'\n');
fprintf(1,['Perform k = %2.0f iterations of AIR Tools sart calling ',...
    'ASTRA for forward and backprojection.\n'],...
    num_iters);
[Xsart_astra,info_astra] = sart(astra_fun_gpu,b,1:num_iters,[],opts);

% ASTRA builtin
fprintf(1,'\n');
fprintf(1,['Perform k = %2.0f iterations of ASTRA''s built-in sirt ',...
    'method, equivalent to AIR Tools sart.\n'],...
    num_iters);
tic_builtin = tic;
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, rot90(B,1));

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% Set up the parameters for a reconstruction algorithm using the GPU
cfg = astra_struct('SIRT_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run the specified number of iterations of the algorithm
astra_mex_algorithm('iterate', alg_id, num_iters);

% Get the result
Xsart_astra_builtin = flipud(astra_mex_data2d('get', rec_id));

% Clean up. 
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
time_builtin = toc(tic_builtin);

% Display all reconstructions
figure
s1 = 2;
s2 = 2;
ca = [0,0.6];

subplot(s1,s2,1)
imagesc(reshape(Xsart_air(:,end),N,N))
axis image off
caxis(ca)
title(sprintf('AIR Tools matrix (%1.2f s)',info_air.timetaken))

subplot(s1,s2,2)
imagesc(reshape(Xsart_air_mf(:,end),N,N))
axis image off
caxis(ca)
title(sprintf('AIR Tools matrix-free (%1.2f s)',info_air_mf.timetaken))

subplot(s1,s2,3)
imagesc(reshape(Xsart_astra(:,end),N,N))
axis image off
caxis(ca)
title(sprintf('AIR Tools calling ASTRA (%1.2f s)',info_astra.timetaken))

subplot(s1,s2,4)
imagesc(reshape(Xsart_astra_builtin,N,N))
axis image off
caxis(ca)
title(sprintf('Built-in ASTRA (%1.2f s)',time_builtin))

colormap gray