function varargout = symkaczmarz(varargin)
%SYMKACZMARZ  Kaczmarz's method with symmetric (top-bottom-top) row sweep
%
%   [X,info] = symkaczmarz(A,b,K)
%   [X,info] = symkaczmarz(A,b,K,x0)
%   [X,info] = symkaczmarz(A,b,K,x0,options)
%
% Implements the symmetric Kaczmarz method: each symkaczmarz sweep 
% consists of a Kaczmarz sweep followed by a Kaczmarz sweep with the 
% equations in reverse order. A symkaczmarz sweep therefore comprises
% approximately double the work of a kaczmarz sweep. To allow fair
% comparison with the other ART methods after a chosen number of 
% iterations, a symkaczmarz sweep is counted as two iterations. To only
% allow full symkaczmarz sweeps, only iterates at even iteration numbers
% can be returned.
%
% Input:
%   A        m times n matrix, or a function that implements matrix-vector
%            multiplication with A and A'; please see explanation below.
%   b        m times 1 vector containing the right-hand side.
%   K        Number of iterations. If K is a scalar, then K is the 
%            maximum number of iterations and only the last iterate is 
%            returned. If K is a vector, then max(K) is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are returned, together with the last iterate.
%            As explained above, only even iteration numbers allowed.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%      relaxpar The relaxation parameter. If relaxpar is a scalar < 2 then
%               it is used in each iteration; default value is 1.
%               Alternatively, relaxpar can be a function with a diminishing
%               parameter, e.g., @(j) 1/sqrt(j), where j counts the total
%               number of row updates.
%      stoprule  Struct containing the following information about the
%                stopping rule:
%                    type = 'none' : (Default) the only stopping rule
%                                    is the maximum number of iterations.
%                           'NCP'  : Normalized Cumulative Perodogram.
%                           'DP'   : Discrepancy Principle.
%                    taudelta   = product of tau and delta, required for DP.
%                    res_dims   = the dimensions that the residual vector
%                                 should be reshaped to, required for NCP.
%                                 E.g. for paralleltomo, res_dims should
%                                 be [p,length(theta)]. For a 1D signal
%                                 res_dims can be a scalar equal to the
%                                 number of elements. 
%                    ncp_smooth = A positive integer specifying the
%                                 filter length in the NCP criterion.
%                                 Default: 2.
%      lbound    Lower bound in box constraint [lbound,ubound]. If scalar,
%                this value is enforced on all elements of x in each 
%                iteration. If vector, it must have same size as x and 
%                then enforces elementwise lower bounds on x. If empty, no
%                bound is enforced. +/-Inf can be used.
%      ubound    Upper bound in box constraint [lbound,ubound]. If scalar,
%                this value is enforced on all elements of x in each 
%                iteration. If vector, it must have same size as x and 
%                then enforces elementwise lower bounds on x. If empty, no
%                bound is enforced. +/-Inf can be used.
%      damp      A parameter P to avoid division by very small row norms
%                by adding P*max_i{||a_i||_2^2} to ||a_i||_2^2.
%      verbose   Nonnegative integer specifying whether progress is printed
%                to screen during iterations. Default=0: no info printed.
%                1: Print in every iteration. Larger than 1: Print every
%                verbose'th iteration and first and last.
%      waitbar   Logical specifying whether a graphical waitbar is shown,
%                default = false.
%
% Output:
%   X        Matrix containing the saved iterations in columns.
%   info     Information struct with 4 fields:
%            stoprule = 0 : stopped by maximum number of iterations
%                       1 : stopped by NCP-rule
%                       2 : stopped by DP-rule
%            finaliter    : no. of iterations in total.
%            relaxpar     : the chosen relaxation parameter.
%            itersaved    : iteration numbers of iterates saved in X.
%            timetaken    : Total time taken by algorithm, in secs.
%
% How to use a function handle for A.
% 1) The user must provide a function myfun that implements matrix-vector
%    multiplication with A and A', with the call
%       y = myfun(v,transp_flag,p1,p2,...)
%    where p1,p2,... are the parameters that define the problem:
%       myfun([],'size',p1,p2,...) returns the size of the matrix,
%       myfun(v,'notransp',p1,p2,...) returns A*v,
%       myfun(w,'transp',p1,p2,...) returns A'*w.
% 2) Before calling symkaczmarz, the user must assign values to the 
%    parameters p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then symkaczmarz is called with this A.
%
% See also: kaczmarz, randkaczmarz, art.

% Reference: Aa. Bjoerck and T. Elfving, Accelerated projection methods for 
% computiong pseudoinverse solutions of systems of linear equations, BIT 
% 19 (1979), pp. 145-163.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

[varargout{1:nargout}] = art('symkaczmarz',varargin{:});