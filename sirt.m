function [X,info,ext_info] = sirt(sirt_method, varargin)
%SIRT General interface for calling SIRT methods.
%
%   [X,info,ext_info] = sirt(sirt_method,A,b,K)
%   [X,info,ext_info] = sirt(sirt_method,A,b,K,x0)
%   [X,info,ext_info] = sirt(sirt_method,A,b,K,x0,options)
%
% Implements a general SIRT method method for the linear system Ax = b:
%
%       x^{k+1} = x^k + relaxpar_k*D*A'*M*(b-A*x^k)
%
% where D and M are diagonal matrices and w_i are weights (default: 
% w_i = 1). All five specific SIRT methods of AIR Tools are special cases
% for specific choices of D and M and computed using this function. In
% addition custom SIRT methods can be specified by other D and M.
%
% Input:
%   sirt_method  Either one of the strings 'landweber', 'cimmino', 'cav',
%                'drop' or 'sart' to specify one of the provided methods.
%                Default is 'sart'.
%                Or a 2-element cell array {Mfun, Dfun} holding two 
%                function handles Mfun and Dfun which implement the
%                multiplication of M and D, respectively. Please see
%                demo_custom for an example.
%   A            m times n matrix, or a function implementing matrix-vector
%                multiplication with A and A'; please see explanation below.
%   b            m times 1 vector containing the right-hand side.
%   K            Number of iterations. If K is a scalar, then K is the 
%                maximum number of iterations and only the last iterate is
%                saved. If K is a vector, then the largest value in K is
%                the maximum number of iterations and only iterates 
%                corresponding to the values in K are saved, together with 
%                the last iterate.
%   x0           n times 1 starting vector. Default: x0 = 0.
%   options      Struct with the following fields:
%      relaxpar  The relaxation parameter. If relaxpar is a scalar then
%                the corresponding value is used in each iteration;
%                default value is 1.9/norm(T*A'*M*A). 
%                If relaxpar is a string, then it refers to a method to 
%                determine relaxpar in each iteration. For this method the
%                following strings can be specified:
%                     'line'    : relaxpar is chosen using line search.
%                     'psi1'    : relaxpar is chosen using the Psi_1-based 
%                                 relaxation method.
%                     'psi1mod' : relaxpar is chosen using the modified 
%                                 Psi_1-based relaxation method.
%                     'psi2'    : relaxpar is chosen using the Psi_2-based
%                                 relaxation method.
%                     'psi2mod' : relaxpar is chosen using the modifed 
%                                 Psi_2-based relaxation method.
%      stoprule  Struct containing the following information about the
%                stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulatice Periodogram.
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
%      s1        Scalar containing largest singular value of sqrt(M)*A.
%      w         m-dimensional weighting vector.
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
%            M            : diagonal of the matrix M = diag(1/||a^i||_S^2).
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
% 2) Before calling sirt, the user must assign values the parameters
%    p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then sirt is called with this A.
%
% For SART using matrix-free it is assumed (not checked) that the
% underlying matrix only has nonnegative elements. If the matrix has one or 
% more negative elements, the result produced by SART is not well-defined. 
% With a sparse matrix, negative elements are allowed and handled properly.
%
% See also: landweber, cimmino, cav, drop, sart.

% Maria Saxild-Hansen, Per Chr. Hansen and Jakob Sauer Jorgensen,
% 2017-03-04 DTU Compute.


% Set default SIRT method to be sart.
if isempty(sirt_method)
    sirt_method = 'sart';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    s1,w,res_dims,rkm1,dk] = check_inputs(varargin{:});

% Extract the Mfun and sfun characterizing each SIRT-type method.
if ischar(sirt_method)
    [Mfun,Dfun] = get_mfun_dfun(sirt_method,varargin{1},m,n,w);
else
    % Possible to pass in custom SIRT method given by 2-element cell array
    % holding function handles to Mfun and Dfun instead of string input.
    Mfun = sirt_method{1};
    Dfun = sirt_method{2};
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialize the values.
k = 0;   % Iteration counter
xk = x0; % Use initial vector.
l = 1;   % Pointing to the next iterate number in K to be saved.

% Do initial check of stopping criteria for x0.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Calculate relaxpar. Special case for SART largest singval is 1.
atma = @(x) Dfun( Afun( Mfun(Afun(x,'notransp')) , 'transp' ));
if strcmpi(sirt_method,'sart')
    s1 = 1;
end
[relaxpar,casel,sigma1tilde] = calc_relaxpar(relaxparinput,s1,kmax,atma,n);

% Store M and D in third output struct if asked for.
if nargout > 2
    ext_info.M = Mfun(ones(m,1));
    ext_info.D = Dfun(ones(n,1));
end

% Main SIRT loop
while ~stop
    
    % Update the iteration number k.
    k = k + 1;

    % Compute the current iteration depending on relaxpar strategy.
    Mrk = Mfun(rk);
    ATMrk = Afun(Mrk,'transp');  %A'*Mrk;  
    if casel == 1
        % SIRT using constant value of relaxpar.
        relaxparcur = relaxpar;        
    elseif casel == 2
        % SIRT using line search.    
        ATMrkS = sum(Dfun(ATMrk.^2));
        relaxparcur = (rk'*Mrk)/ATMrkS;
    elseif casel == 3
        % SIRT using psi1 or psi2.
        relaxparcur = relaxpar(k);
    end % end the different cases of relaxpar strategies.
    
    % The update step with current relaxpar
    xk = xk + relaxparcur*(Dfun(ATMrk));
    
    % Enforce any lower and upper bounds (scalars or xk-sized vectors)
    if ~isempty(lbound)
        xk = max(xk,lbound);
    end
    if ~isempty(ubound)
        xk = min(xk,ubound);
    end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, relaxparcur, taudelta, k, kmax, rkm1, dk, res_dims);
        
    % If the current iteration is requested saved.
    if k == K(l) || stop
        X(:,l) = xk;
        l = l + 1;
    end
end

% Return only the saved iterations: Only to "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);

% Save further to info:

% Largest singular value determined
info.s1 = sigma1tilde;

% List of iterates saved: all in K smaller than the final, and the final.
info.itersaved = [K(K<info.finaliter), info.finaliter];
