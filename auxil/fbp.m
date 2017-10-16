function x = fbp(A,b,theta,filter)
%FBP  Filtered Back Projection using the system matrix
%
% x = fbp(A,b,theta,filter)
%
% Input: A        System matrix,
%        b        Right-hand side (vectorized sinogram),
%        theta    Vector of projection angles,
%        filter   The following filters are available:
%                  'Ram-Lak' (default), 'Shepp-Logan', 'cosine',
%                  'Hamming', 'Hann', 'none' (gives pure backprojection).
%
% Output: x       The solution vector (vectorized reconstruction)
%
% This function computes a filtered back projection solution to a parallel-
% beam tomography problem, such as arising from the AIR Tools II test
% problem paralleltomo.  It is somewhat similar to Matlab's iradon, but it
% conforms with AIR Tools II in that the backprojection is performed by
% multiplication with the transpose of A.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Input check and set default filter.
if nargin < 3, error('Input A, b, theta must me specified'); end
if nargin < 4, filter = 'ram-lak'; end
filter = lower(filter);  % Turn all characters to lower case.

% Pure back projection.
if strcmp(filter,'none'), x = A'*b; return, end

% Reshape right-hand side b to a sinogram.
m = size(A,1);        % Number of rows in A.
na = length(theta);   % Number of angles.
nr = m/na;            % Number of rays.
b = reshape(b,nr,na); % Reshape sinogram, overwrite to save memory.

% Set ramp filter and frequencies.
order = max(64,2^nextpow2(2*nr));
n = (0:(order/2))';
f = n/max(n);
f(1) = 0.4/order;
f(end) = 1 - f(1);
w = 2*pi*n/order;

% Introduce low-pass filter and symmetrice.
switch filter
    case 'ram-lak'
        % No low-pass filtering.
    case 'shepp-logan'
        f(2:end) = f(2:end).*sinc(w(2:end)/(2*pi));
    case 'cosine'
        f(2:end) = f(2:end).*cos(w(2:end)/2);
    case 'hamming'
        f(2:end) = f(2:end).*(.54 + .46*cos(w(2:end)));
    case 'hann'
        f(2:end) = f(2:end).*(1+cos(w(2:end))) / 2;
    otherwise
        error('Iinvalid filter')
end
f = [f ; f(end-1:-1:2)];

% Zero paddding of sinogram, and filtering of each column.
b(length(f),1) = 0;
b = fft(b);
for i = 1:na
    b(:,i) = b(:,i).*f; % frequency domain filtering
end
b = real(ifft(b));
b(nr+1:end,:) = [];

% Backprojection.
x = A'*b(:)*pi/(2*length(theta));