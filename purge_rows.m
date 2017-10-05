function [A,b] = purge_rows(A,b,Nthr)
%PURGE_ROWS  Remove zero or very sparse rows of A and corresp. entries in b
%
%   [A,b] = purge_rows(A,b)
%   [A,b] = purge_rows(A,b,Nthr)
%
% Identifies zero rows of the coefficient matrix A and removes them.
% If a right-hand side b is present, the corresponding elements of b are
% also removed (b can be a matrix of several right hand sides). A must be
% matrix, i.e., it cannot be a function handle.
%
% If a positive Nthr is given as the third argument, then all rows for
% which the number of nonzeros is less than or equal to Nthr are removed.
%
% Use this function to 'clean up' a discretized tomography problem.
% Zero rows do not contribute to the reconstruction.
% Rows with few nonzero elements correspond to pixels near the corners of
% the image, whose reconstructions are highly sensitive to noise.
%
% See also: demo_custom_sirt.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

if nargin<3, Nthr = 0; end

s = sum(A~=0,2);
I = find(s>Nthr);
A = A(I,:);

if nargin>1 && ~isempty(b), b = b(I,:); end