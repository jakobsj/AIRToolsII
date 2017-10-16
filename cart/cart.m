function [X,info] = cart(cart_method, varargin)
%CART  General interface for the CART methods
%
%   [X,info] = cart(cart_method,A,b,K)
%   [X,info] = cart(cart_method,A,b,K,x0)
%   [X,info] = cart(cart_method,A,b,K,x0,options)
%
% Implements a general CART method for the system Ax = b:
%
%       x_j^{k+1} = x_j^k + relaxpar*a_j'(b - A x^k)/(||a_j||_2^2)
%
% where a_j is the j-th column of A, and the order for j can be chosen in 
% different ways. The provided method columnaction uses the order from
% first to last column, while other orders can be specified by the user. 
%
% This method incorporates the "flagging" idea: a component of the
% solution x is "flagged" if its update is smaller than a threshold THR
% times max(abs(x)) is not updated again until it is "unflagged" which
% happens after a random number of iterations chosen as rand*Nunflag.
%
% Input:
%   cart_method  Either the string 'columnaction' to specify the provided 
%                method; default is 'columnaction'.
%                Or a column index vector of length n (for matrix A m-by-n)
%                with a desired fixed order in which to step through all 
%                columns of A. Please see demo_custom for an example.
%   A            m times n matrix, or a function that implements matrix-
%                vector multiplication with A and A'; see explanation below.
%   b            m times 1 vector containing the right-hand side.
%   K            Number of iterations. If K is a scalar, then K is the 
%                maximum number of iterations and only the last iterate is 
%                returned. If K is a vector, then max(K) is the maximum
%                number of iterations and only iterates corresponding to the
%                values in K are returned, together with the last iterate.
%   x0           n times 1 starting vector. Default: x0 = 0.
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
%                    ncp_smooth = An positive integer specifying the
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
%      THR       A component is "flagged" if its update is smaller than
%                THR*max(abs(x)) where x is the previous iteration vector.
%                Default = 1e-4.
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
% 2) Before calling cart, the user must assign values to the 
%    parameters p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then cart is called with this A.
%
% See also: art, columnaction.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.
% With contribution from Jacob Frosig and Nicolai Riis.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Measure total time taken.
t_total = tic;

% Set default CART method to be columnaction.
if isempty(cart_method)
    cart_method = 'columnaction';
end

% Parse inputs.
[Afun,b,~,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    ~,res_dims,rkm1,dk,do_waitbar,verbose,damp,THR,Kbegin,Nunflag] = ...
    check_inputs(varargin{:});

% Faster to access rows of matrix directly if available.
A = varargin{1};

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialization before iterations.
k = 0;   % Iteration counter
xk = x0; % Use initial vector.
l = 1;   % Pointing to the next iterate number in K to be saved.

% Do initial check of stopping criteria - probably relaxpar should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Compute the relaxation parameter to be used throughout iterations.
relaxpar = calc_relaxpar(relaxparinput);

% Calculate the norm of each column in A.
normAj = zeros(1,n);
if isnumeric(A)
    B = 200;  % Block size; adjust if necessary.
    for I = 1:ceil(n/B)
        J = 1+(I-1)*B : min(n,I*B);
        normAj(J) = sum(A(:,J).*A(:,J),1);
    end
elseif isa(A,'function_handle')
    for j = 1:n
        e = zeros(n,1); e(j) = 1;
        v = A(e,'notransp');
        normAj(j) = norm(v)^2;
    end
else
    for j = 1:n
        e = zeros(n,1); e(j) = 1;
        v = A*e;
        normAj(j) = norm(v)^2;
    end
end

% Depending on CART method, set the column order.
if ischar(cart_method)
    switch lower(cart_method)
        case 'columnaction'
            % Only loop over nonzero columns.
            J = find(normAj>0);
        otherwise
            error('Unknown CART method specified')
    end
else
    % Custom CART method specified by user giving the column order J as
    % first input instead of the string with a particular CART method name.
    J = cart_method(:)';
    J = J(normAj(J)>0);
end

% Apply damping.
normAj = normAj + damp*max(normAj);

% Set additional CART "flagging" parameters.
F = true(n,1);       % Vector of logical "flags."
Nflag = zeros(n,1);  % Counts "flagged" iterations for each component.

% For deciding how to apply constraints. Can be nan, scalar or vector.
is_lbound_empty  = isempty(lbound);
is_lbound_scalar = isscalar(lbound);
is_ubound_empty  = isempty(ubound);
is_ubound_scalar = isscalar(ubound);

% Initalize waitbar if selected.
if do_waitbar
    h_waitbar = waitbar(0);
end

% Main CART loop.
while ~stop
            
    % Update timer for current iteration
    t_iter = tic;
    
    % Update the iteration number k.
    k = k + 1;
    
    % Update waitbar if selected
    if do_waitbar
        waitbar(k/kmax,h_waitbar,...
            sprintf('Running iteration %d of %d...',k, kmax))
    end
    
    % The columnwise sweep.
    for j = J
        mm = max(abs(xk)); % Max abs element of previous iteration.
        if F(j)
            
            % Get the j'th column of A
            if isnumeric(A)
                aj = A(:,j);
            elseif isa(A,'function_handle')
                e = zeros(n,1); e(j) = 1;
                aj = A(e,'notransp');
            else
                e = zeros(n,1); e(j) = 1;
                aj = A*e;
            end
            
            % The update.
            delta = aj'*rk/normAj(j);
            od = relaxpar*delta;           
        
            % Correction for constraints, first get current value for reuse
            xkj = xk(j);
            
            % Apply, only if not NaN, and if od would take xkj below lbj or
            % above ubj, correction to od, so new xk(j) will be set to
            % either lbj or ubj.
            if ~is_lbound_empty
                % Handle that lbound can be either scalar or vector: If
                % scalar extract first element, if non-scalar extract
                % element j of lbound.
                lbj = lbound( (~is_lbound_scalar)*j + is_lbound_scalar );
                if od < lbj - xkj
                    od = lbj - xkj;
                end
            end
            % Same for ubound.
            if ~is_ubound_empty
                % Same for ubound.
                ubj = ubound( (~is_ubound_scalar)*j + is_ubound_scalar );
                if od > ubj - xkj
                    od = ubj - xkj;
                end
            end
            
            % Apply update.
            xk(j) = xkj + od;
            
            % "Flag" if needed.
            if k >= Kbegin && abs(od) < THR*mm;
                F(j) = false;
                Nflag(j) = 1;  % PCH
            end
            
            % Update residual.
            rk = rk - od*aj;
        else
            if Nflag(j) < randi(Nunflag)
                Nflag(j) = Nflag(j) + 1;
            else
                F(j) = true;
            end
        end
    end
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims);
    
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
            k,kmax,norm(rk),relaxpar,toc(t_iter),toc(t_total));
    end
end

% Close waitbar if selected.
if do_waitbar
    close(h_waitbar);
end

% Return only the saved iterations: use "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);

% List of iterates saved: all in K smaller than the final, and the final.
K = K(:);  info.itersaved = [K(K<info.finaliter); info.finaliter];

% Save time total time taken
info.timetaken = toc(t_total);
