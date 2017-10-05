function varargout = cav(varargin)
%CAV  Component Averaging (CAV) method
%
%   [X,info,ext_info] = cav(A,b,K)
%   [X,info,ext_info] = cav(A,b,K,x0)
%   [X,info,ext_info] = cav(A,b,K,x0,options)
%
% Implements the CAV method for the linear system Ax = b:
%
%       x^{k+1} = x^k + relaxpar_k*A'*M*(b-A*x^k)
%
% where M = diag(1/||a_i||_S^2, S = diag(s_j), and s_j denotes the number
% of nonzero elements in column j.
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
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%      relaxpar  The relaxation parameter. If relaxpar is a scalar then
%                the corresponding value is used in each iteration;
%                default value is 1.9/norm(A'*M*A). 
%                If relaxpar is a string, then it refers to a method to 
%                determine relaxpar in each iteration. For this method the
%                following strings can be specified:
%                    'line'    : relaxpar is chosen using line search.
%                    'psi1'    : relaxpar is chosen using the Psi_1-based 
%                                relaxation method.
%                    'psi1mod' : relaxpar is chosen using the modified 
%                                Psi_1-based relaxation method.
%                    'psi2'    : relaxpar is chosen using the Psi_2-based
%                                relaxation method.
%                    'psi2mod' : relaxpar is chosen using the modifed 
%                                Psi_2-based relaxation method.
%      stoprule  Struct containing the following information about the
%                stopping rule:
%                    type = 'none' : (Default) the only stopping rule
%                                    is the maximum number of iterations.
%                           'NCP'  : Normalized Cumulatice Periodogram.
%                           'DP'   : Discrepancy Principle.
%                           'ME'   : Monotone Error rule.
%                    taudelta   = product of tau and delta, required for
%                                 DP and ME.
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
%      rho       Scalar containing spectral radius of the iteration matrix.
%      verbose   Nonnegative integer specifying whether progress is printed
%                to screen during iterations. Default=0: no info printed.
%                1: Print in every iteration. Larger than 1: Print every
%                verbose'th iteration and first and last.
%      waitbar   Logical specifying whether a graphical waitbar is shown,
%                default = false.
%
% Output:
%   X        Matrix containing the saved iterations in columns.
%   info     Information struct with 5 fields:
%            stoprule = 0 : stopped by maximum number of iterations
%                       1 : stopped by NCP-rule
%                       2 : stopped by DP-rule
%                       3 : stopped by ME-rule.
%            finaliter    : no. of iterations in total.
%            relaxpar     : the chosen relaxation parameter.
%            rho          : the computed spectral radius.
%            itersaved    : iteration numbers of iterates saved in X.
%            timetaken    : Total time taken by algorithm, in secs.
%   ext_info Extra information struct with 2 fields:
%            M            : diagonal of the matrix M = diag(1/||a_i||_S^2).
%            D            : diagonal of the matrix D (all ones for cav)
%
% How to use a function handle for A.
% 1) The user must provide a function myfun that implements matrix-vector
%    multiplication with A and A', with the call
%       y = myfun(v,transp_flag,p1,p2,...)
%    where p1,p2,... are the parameters that define the problem:
%       myfun([],0,p1,p2,...) returns the size of the matrix,
%       myfun(v,'notransp',p1,p2,...) returns A*v,
%       myfun(w,'transp',p1,p2,...) returns A'*w.
% 2) Before calling cav, the user must assign values the parameters
%    p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then cav is called with this A.
%
% See also: cimmino, drop, landweber, sart, sirt.

% Reference: Y. Censor, D. Gordon, and R. Gordon, Component averaging: An 
% efficient iterative parallel algorithm for large sparse unstructured 
% problems, Parallel Computing, 27 (2001), pp. 777-808.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

[varargout{1:nargout}] = sirt('cav',varargin{:});
