%DEMO_SHOW_TOMO  Illustrates the use of the show_tomo function
%
% This script illustrates each of the 6 supplied tomographic test problems.
% First, while generating the test problem, the ray paths used in each
% problem are plotted. Second, after generation of the test problem is
% complete, the rows of the generated matrix are displayed as an image to
% illustrate which entries in the system matrix arise in each geometry.
% For the parallel-beam geometry we also illustration the corresponding use
% of the matrix-free function handle representation of the test problem.
%
% See also: paralleltomo, fancurvedtomo, fanlineartomo, show_tomo,
% sphericaltomo, seismictomo, seismicwavetomo.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

clear, clc
fprintf(1,'Starting demo_show_tomo:\n\n');

% Set the parameters for the test problem.
N = 32;           % The discretization points.
timedelay = 0.1;  % Time between display of rays.

%% Paralleltomo, matrix and matrixfree - parallel-beam geometry

fprintf(1,'paralleltomo...\n\n');

theta = 0:30:150; % No. of used angles.
p = 32;           % No. of parallel rays.

% Create matrix version of test problem and display rays 
A = paralleltomo(N,theta,p,[],timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)

% Create the function handle version, this time omitting display of rays.
AF = paralleltomo(N,theta,p,[],0,0);

% Illustrate the rows of the function handle A
figure
show_tomo(AF,[],timedelay)


%% Fancurvedtomo - fan-beam geometry with curved detector

fprintf(1,'fancurvedtomo...\n\n');

theta = 0:60:300; % No. of used angles.
p = 32;           % No. of parallel rays.
R = 2;            % Source to center distance, in units of N.
d = 28;           % Angular range of rays, in degrees

% Create matrix version of test problem and display rays 
A = fancurvedtomo(N,theta,p,R,d,timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)


%% Fanlineartomo - fan-beam geometry with linear detector

fprintf(1,'fanlineartomo...\n\n');

theta = 0:60:300; % No. of used angles.
p = 32;           % No. of parallel rays.
R = 2;            % Source to center distance, in units of N.
sd = 3;           % Source to detector distance, in units of N.
dw = sd/R;        % Detector size

% Create matrix version of test problem and display rays 
A = fanlineartomo(N,theta,p,R,dw,sd,timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)


%% Sphericaltomo - 2D spherical Radon tomography test problem

fprintf(1,'sphericaltomo...\n\n');

theta = 0:60:300; % No. of used angles.
numCircles = 32;  % No. of circles at each angular position.

% Create matrix version of test problem and display rays 
A = sphericaltomo(N,theta,numCircles,timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)


%% Seismictomo - 2D seismic travel-time tomography test problem

fprintf(1,'seismictomo...\n\n');

s = 6;   % The number of sources in the right side of the domain.
p = 10;  % The number of receivers (seismographs) equally spaced on the
         % surface and on the left side of the domain (default p = 2*N).

% Create matrix version of test problem and display rays 
A = seismictomo(N,s,p,timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)


%% Seismicwavetomo - Seismic tomography problem w/o ray assumption

fprintf(1,'seismicwavetomo...\n\n');

s = 6;      % The number of sources in the right side of the domain.
p = 10;     % The number of receivers (seismographs) equally spaced on the
            % surface and on the left side of the domain (default p = 2*N).
omega = 10; % Dominant frequency of the propagating wave (default = 10).

% Create matrix version of test problem and display rays 
A = seismicwavetomo(N,s,p,omega,timedelay);

% Illustrate the rows of the matrix A.
figure
show_tomo(A,[],timedelay)