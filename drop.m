function [X,info,ext_info] = drop(varargin)
%DROP Diagonally Relaxed Orthogonal Projections (DROP) method
%
%   [X,info,ext_info] = drop(A,b,K)
%   [X,info,ext_info] = drop(A,b,K,x0)
%   [X,info,ext_info] = drop(A,b,K,x0,options)
%
% Implements the DROP method for the linear system Ax = b:
%
%       x^{k+1} = x^k + relaxpar_k*D*A'*M*(b-A*x^k)
%
% where D = diag(1/s_j), s_j is the number of nonzero elements in the j'th
% column of A, M = diag(w_i/||a^i||_2^2), and w_i are weights (default w_i = 1).
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
%       relaxpar  The relaxation parameter. If relaxpar is a scalar then
%                 the corresponding value is used in each iteration;
%                 default value is 1.9/norm(T*A'*M*A). 
%                 If relaxpar is a string, then it refers to a method to 
%                 determine relaxpar in each iteration. For this method the
%                 following strings can be specified:
%                     'line'    : relaxpar is chosen using line search.    
%                     'psi1'    : relaxpar is chosen using the Psi_1-based 
%                                 relaxation method.
%                     'psi1mod' : relaxpar is chosen using the modfied 
%                                 Psi_1-based relaxation method.
%                     'psi2'    : relaxpar is chosen using the Psi_2-based
%                                 relaxation method.
%                     'psi2mod' : relaxpar is chosen using the modfied 
%                                 Psi_2-based relaxation method.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Periodogram.
%                            'DP'   : Discrepancy Principle.
%                            'ME'   : Monotone Error rule.
%                     taudelta   = product of tau and delta, only needed
%                                  for DP and ME.
%                     res_dims   = the dimensions that the residual vector
%                                  should be reshaped to, required for NCP.
%                                  E.g. for paralleltomo, res_dims should
%                                  b e [p,length(theta)]. For a 1D signal
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
%       s1        Scalar containing largest singular value of sqrt(M)*A.
%       w         m-dimensional weighting vector.
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
%            s1           : the computed largest singular value.
%            itersaved    : iteration numbers of iterates saved in X.
%   ext_info Extra information struct with 2 fields:
%            M            : diagonal of the matrix M = diag(w_i/||a^i||^2_2)
%            D            : diagonal of the matrix D = diag(1/s_j).
%
% How to use a function handle for A.
% 1) The user must provide a function myfun that implements matrix-vector
%    multiplication with A and A', with the call
%       y = myfun(v,transp_flag,p1,p2,...)
%    where p1,p2,... are the parameters that define the problem:
%       myfun([],0,p1,p2,...) returns the size of the matrix,
%       myfun(v,'notransp',p1,p2,...) returns A*v,
%       myfun(w,'transp',p1,p2,...) returns A'*w.
% 2) Before calling landweber, the user must assign values the parameters
%    p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then drop is called with this A.
%
% See also: landweber, cimmino, cav, sart.

% Maria Saxild-Hansen, Per Chr. Hansen and Jakob Sauer Jorgensen,
% November 8, 2015, DTU Compute.

% Reference: Y. Censor, T. Elfving, G. T. Herman, and T. Nikazad, On
% diagonally relaxed orthogonal projection methods, SIAM J. Sci. Comp.,
% 30 (2007), pp. 473-504.

[X,info,ext_info] = sirt('drop',varargin{:});