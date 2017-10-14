function varargout = columnaction(varargin)
%COLUMNACTION  Column-action method with cyclic column sweeps
%
%   [X,info] = columnaction(A,b,K)
%   [X,info] = columnaction(A,b,K,x0)
%   [X,info] = columnaction(A,b,K,x0,options)
%
% Implements the column-action method for the system Ax = b:
%
%       x_j^{k+1} = x_j^k + relaxpar*a_j'(b - A x^k)/(||a_j||_2^2)
%
% where a_j is the j-th column of A, and j = (k mod n) + 1.
%
% This method is identical to the coordinate descent optimization method.
% Our implementation incorporates the "flagging" idea: a component of the
% solution x is "flagged" if its update is smaller than a threshold THR
% times max(abs(x)); it is not updated again until it is "unflagged" which
% happens after a random number of iterations chosen as rand*Nunflag.
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
%      relaxpar  The relaxation parameter. For this method relaxpar must
%                be a scalar; default value is 0.25.
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
%                by adding P*max_i{||a_j||_2^2} to ||a_j||_2^2.
%      THR       A component is "flagged" if its update is smaller than
%                THR*max(abs(x)) where x is the previous iteration vector.
%                Default THR = 1e-4.
%      Kbegin    Perform Kbegin iterations before allowing "flagging".
%                Default = 10.
%      Nunflag   Upper bound on number of iterations to wait before
%                unflagging a component. The number of iterations to wait
%                is chosen as randi(Nunflag). Default = round(max(K/4).
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
% 2) Before calling columnaction, the user must assign values to the 
%    parameters p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then columnaction is called with this A.
%
% See also: cart, art, kaczmarz.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

[varargout{1:nargout}] = cart('columnaction',varargin{:});
