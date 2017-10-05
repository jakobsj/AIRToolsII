function [X,info,ext_info] = sirt(sirt_method, varargin)
%SIRT  General interface for all SIRT methods
%
%   [X,info,ext_info] = sirt(sirt_method,A,b,K)
%   [X,info,ext_info] = sirt(sirt_method,A,b,K,x0)
%   [X,info,ext_info] = sirt(sirt_method,A,b,K,x0,options)
%
% Implements a general SIRT method method for the linear system Ax = b:
%
%       x^{k+1} = x^k + relaxpar_k*D*A'*M*(b-A*x^k)
%
% where D and M are matrices, typically diagonal, but non-diagonal also
% allowed. All five specific SIRT methods of AIR Tools are special cases
% for specific choices of D and M and computed using this function. In
% addition custom SIRT methods can be specified by other D and M.
%
% Input:
%   sirt_method  Either one of the strings 'landweber', 'cimmino', 'cav',
%                'drop' or 'sart' to specify one of the provided methods;
%                default is 'sart'.
%                Or a struct with fields M and D, each holding either a
%                matrix or a vector representing the diagonal elements of a
%                diagonal matrix of appropriate size. In this way, new and 
%                completely general SIRT methods can be specified by the
%                choice of the M and D matrices. If either is omitted, an
%                identity matrix is used be default. Please see
%                demo_custom_all for an example.
%   A            m times n matrix, or a function that implements matrix-
%                vector multiplication with A and A'; see explanation below.
%   b            m times 1 vector containing the right-hand side.
%   K            Number of iterations. If K is a scalar, then K is the 
%                maximum number of iterations and only the last iterate is 
%                returned. If K is a vector, then max(K) is the maximum
%                number of iterations and only iterates corresponding to the
%                values in K are returned, together with the last iterate.
%   x0           n times 1 starting vector. Default: x0 = 0.
%   options      Struct with the following fields:
%      relaxpar  The relaxation parameter. If relaxpar is a scalar then
%                it is used in each iteration; default = 1.9/rho(T*A'*M*A). 
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
%                     taudelta   = product of tau and delta, required for
%                                  DP and ME.
%                     res_dims   = the dimensions that the residual vector
%                                  should be reshaped to, required for NCP.
%                                  E.g. for paralleltomo, res_dims should
%                                  be [p,length(theta)]. For a 1D signal
%                                  res_dims can be a scalar equal to the
%                                  number of elements. 
%                     ncp_smooth = A positive integer specifying the
%                                  filter length in the NCP criterion.
%                                  Default: 2.
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
%            M            : Diagonal of matrix M (if diagonal) or matrix.
%            D            : Diagonal of matrix D (if diagonal) or matrix.
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
% For SART using a function handle, it is assumed (not checked) that the
% underlying matrix only has nonnegative elements. If the matrix has one or 
% more negative elements, the result produced by SART is not well-defined. 
% With a sparse matrix, negative elements are allowed and handled properly.
%
% See also: landweber, cimmino, cav, drop, sart.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Measure total time taken.
t_total = tic;

% Set default SIRT method to be sart.
if isempty(sirt_method)
    sirt_method = 'sart';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    rho,res_dims,rkm1,dk,do_waitbar,verbose] = check_inputs(varargin{:});

% Extract the Mfun and Dfun characterizing each SIRT-type method.
[Mfun,Dfun,Mflag,Dflag] = get_mfun_dfun(sirt_method,varargin{1},m,n);

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
    rho = 1;
end
[relaxpar,casel,rho] = calc_relaxpar(relaxparinput,rho,kmax,atma,n);

% Store Mfun and Dfun in third output struct if asked for.
if nargout > 2
    if Mflag == 0
        % M is identity matrix implemented matrix-free. Create vector
        % holding 1-values from diagonal.
        ext_info.M = ones(m,1);
    else
        % M is a vector or matrix and can be extracted from function handle
        % by scalar multiplication by 1.
        ext_info.M = Mfun(1);
    end
    if Dflag == 0
        % D is identity matrix implemented matrix-free. Create vector
        % holding 1-values from diagonal.
        ext_info.D = ones(n,1);
    else
        % D is a vector or matrix and can be extracted from function handle
        % by scalar multiplication by 1.
        ext_info.D = Dfun(1);
    end
end

% Initalize waitbar if selected.
if do_waitbar
    h_waitbar = waitbar(0);
end

% Main SIRT loop
while ~stop
    
    % Update timer for current iteration.
    t_iter = tic;
    
    % Update the iteration number k.
    k = k + 1;
    
    % Update waitbar if selected.
    if do_waitbar
        waitbar(k/kmax,h_waitbar,...
            sprintf('Running iteration %d of %d...',k, kmax))
    end

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
    
    % Enforce any lower and upper bounds (scalars or xk-sized vectors).
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
    
    % Print info, if selected. If verbose is 0, do not print. If 1, print
    % at every iteration. If larger than 1, print first, every verbose'th
    % iteration and final iteration when stopping rule is met, including 
    % if reaching maximum iterations.
    if verbose == 1 || (verbose>1 && (k==1 || stop || mod(k,verbose)==1))
        fprintf(['%*d/%d:  resnorm=%1.2e  relaxpar=%1.2e  ',...
            'time(iter/total)=%1.2e/%1.2e secs.\n'], ceil(log10(kmax)), ...
            k,kmax,norm(rk),relaxparcur,toc(t_iter),toc(t_total));
    end
end

% Close waitbar if selected.
if do_waitbar
    close(h_waitbar);
end

% Return only the saved iterations: use "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);

% Save further to info:

% Largest singular value determined
info.rho = rho;

% List of iterates saved: all in K smaller than the final, and the final.
K = K(:); info.itersaved = [K(K<info.finaliter); info.finaliter];

% Save time total time taken
info.timetaken = toc(t_total);
