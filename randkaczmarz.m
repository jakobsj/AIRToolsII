function [X,info] = randkaczmarz(varargin)
%RANDKACZMARZ Randomized Kaczmarz method
%
%   [X,info] = randkaczmarz(A,b,K)
%   [X,info] = randkaczmarz(A,b,K,x0)
%   [X,info] = randkaczmarz(A,b,K,x0,options)
%
% Implements the randomized Kaczmarz's method for the system Ax = b:
%
%       x^{k+1} = x^k + relaxpar*(b_i - a^i*x^k)/(||a^i||_2^2)a^i
%
% where row i is chosen randomly with probability proportional to ||a^i||_2^2.
% One iteration consists of m such steps.
%
% Input:
%   A        m times n matrix, or a function that implements matrix-vector
%            multiplication with A and A'; please see explanation below.
%   b        m times 1 vector containing the right-hand side.
%   K        Number of iterations. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are saved, together with the last iterate.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       relaxpar  The relaxation parameter. For this method relaxpar must
%                 be a scalar < 2; default value is 1.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Perodogram.
%                            'DP'   : Discrepancy Principle.
%                     taudelta   = product of tau and delta, only needed
%                                  for DP.
%                     res_dims   = the dimensions that the residual vector
%                                  should be reshaped to, required for NCP.
%                                  E.g. for paralleltomo, res_dims should
%                                  be [p,length(theta)]. For a 1D signal
%                                  res_dims can be a scalar equal to the
%                                  number of elements. 
%                     ncp_smooth = An positive integer specifying number of
%                                  iterations to filter/average NCP
%                                  criterion over. Default: 4.
%       lbound    Lower bound in box constraint [lbound,ubound]. If scalar,
%                 this value is enforced on all elements of x in each 
%                 iteration. If vector, it must have same size as x and 
%                 then enforces elementwise lower bounds on x. If empty, no
%                 bound is enforced. +/-Inf can be used.
%       ubound    Upper bound in box constraint [lbound,ubound]. If scalar,
%                 this value is enforced on all elements of x in each 
%                 iteration. If vector, it must have same size as x and 
%                 then enforces elementwise lower bounds on x. If empty, no
%                 bound is enforced. +/-Inf can be used.
%       damping   A parameter P to avoid division by very small row norms
%                 by adding P*max_i{||a^i||_2^2} to ||a^i||_2^2.
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
%
% How to use a function handle for A.
% 1) The user must provide a function myfun that implements matrix-vector
%    multiplication with A and A', with the call
%       y = myfun(v,transp_flag,p1,p2,...)
%    where p1,p2,... are the parameters that define the problem:
%       myfun([],'size',p1,p2,...) returns the size of the matrix,
%       myfun(v,'notransp',p1,p2,...) returns A*v,
%       myfun(w,'transp',p1,p2,...) returns A'*w.
% 2) Before calling randkaczmarz, the user must assign values to the 
%    parameters p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then randkaczmarz is called with this A.
%
% See also: kaczmarz, symkaczmarz, art.

% Maria Saxild-Hansen and Per Chr. Hansen, Nov. 8, 2015, DTU Compute.

% Reference: T. Strohmer and R. Vershynin, A randomized Kaczmarz algorithm
% for linear systems with exponential convergence, J. Fourier Analysis and
% Applications, 15 (2009), pp. 262-278.

[X,info] = art('randkaczmarz',varargin{:});