function [X,info] = art(art_method, varargin)
%ART  General interface for all Kaczmarz/ART methods
%
%   [X,info] = art(art_method,A,b,K)
%   [X,info] = art(art_method,A,b,K,x0)
%   [X,info] = art(art_method,A,b,K,x0,options)
%
% Implements a general ART method for the linear system Ax = b:
%       
%       x^{k+1} = x^k + relaxpar*(b_i - a_i'*x^k)/(||a_i||_2^2)*a_i
%
% where a_i' is the i-th row of A, and the order for i is chosen in 
% different ways by the provided methods (see kaczmarz, symkaczmarz, 
% randkaczmarz for details) or a custom order specified by the user.
% One iteration consists of m such steps (m-1 for symkaczmarz).
%
% Input:
%   art_method  Either one of the strings 'kaczmarz', 'symkaczmarz', or 
%               'randkaczmarz' to specify one of the provided methods;
%               default is 'kaczmarz'.
%               Or a row index vector of length m (for matrix A m-by-n)
%               with a desired fixed order in which to step through all 
%               rows of A. Please see demo_custom for an example.
%   A           m times n matrix, or a function that implements matrix-
%               vector multiplication with A and A'; see explanation below.
%   b           m times 1 vector containing the right-hand side.
%   K           Number of iterations. If K is a scalar, then K is the 
%               maximum number of iterations and only the last iterate is 
%               returned. If K is a vector, then max(K) is the maximum
%               number of iterations and only iterates corresponding to the
%               values in K are returned, together with the last iterate.
%   x0          n times 1 starting vector. Default: x0 = 0.
%   options     Struct with the following fields:
%      relaxpar The relaxation parameter. If relaxpar is a scalar < 2 then
%               it is used in each iteration; default value is 1.
%               Alternatively, relaxpar can be a function with a diminishing
%               parameter, e.g., @(j) 1/sqrt(j), where j counts the total
%               number of row updates.
%      stoprule Struct containing the following information about the
%               stopping rule:
%                   type = 'none' : (Default) the only stopping rule
%                                   is the maximum number of iterations.
%                          'NCP'  : Normalized Cumulative Perodogram.
%                          'DP'   : Discrepancy Principle.
%                   taudelta   = product of tau and delta, required for DP.
%                   res_dims   = the dimensions that the residual vector
%                                should be reshaped to, required for NCP.
%                                E.g., for paralleltomo res_dims should
%                                be [p,length(theta)]. For a 1D signal,
%                                res_dims can be a scalar equal to the
%                                number of elements. 
%                   ncp_smooth = An positive integer specifying the filter
%                                length in the NCP criterion; default is 2.
%      lbound   Lower bound in box constraint [lbound,ubound]. If scalar,
%               this value is enforced on all elements of x in each 
%               iteration. If vector, it must have same size as x and 
%               then enforces elementwise lower bounds on x. If empty, no
%               bound is enforced. +/-Inf can be used.
%      ubound   Upper bound in box constraint [lbound,ubound]. If scalar,
%               this value is enforced on all elements of x in each 
%               iteration. If vector, it must have same size as x and 
%               then enforces elementwise lower bounds on x. If empty, no
%               bound is enforced. +/-Inf can be used.
%      damp     A parameter damp to avoid division by very small row norms
%               by adding damp*max_i{||a_i||_2^2} to ||a_i||_2^2.
%      verbose  Nonnegative integer specifying whether progress is printed
%               to screen during iterations. Default=0: no info printed.
%               1: Print in every iteration. Larger than 1: Print every
%               verbose'th iteration and first and last.
%      waitbar  Logical specifying whether a graphical waitbar is shown,
%               default = false.
%
% Output:
%   X        Matrix containing the saved iterations as the columns.
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
% 2) Before calling art, the user must assign values the parameters
%    p1,p2,... and define an new function handle A in this way:
%       A = @(v,transp_flag) myfun(v,transp_flag,p1,p2,...);
% 3) Then art is called with this A.
%
% See also: cart, kaczmarz, randkaczmarz, symkaczmarz.

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

% Set default ART method to be kaczmarz.
if isempty(art_method)
    art_method = 'kaczmarz';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
     ~,res_dims,rkm1,dk,do_waitbar,verbose,...
     damp] = check_inputs(varargin{:});

% Special check for symkaczmarz: number of iterations must be even.
if ischar(art_method) && strncmpi(art_method,'sym',3)
    if any(mod(K,2))
        error('For symkaczmarz only even iteration numbers can be requested.');
    else
        % Since in symkaczmarz a sweep is top -> bottom -> top, it is twice
        % as expensive as other methods, and one sweep is therefore counted
        % as two iterations. For code reasons, we do this by running half
        % the number of iterations, then later doubling iteration numbers
        % again.
        K = K/2;
        kmax = kmax/2;
    end    
end

% If A is matrix, it is more efficient for ART to work with transposed A.
A = varargin{1};
if ~isa(A,'function_handle')
    A = A';
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialization before iterations.
k = 0;   % Iteration counter.
xk = x0; % Use initial vector.
l = 1;   % Pointing to the next iterate number in K to be saved.

% Initial check of stopping criteria.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Compute the relaxation parameter to be used throughout iterations.
relaxpar = calc_relaxpar(relaxparinput);

% Calculate the norm of each row in A. Remember that A is transposed.
normAi = zeros(1,m);
if isnumeric(A)
    B = 200;  % Block size; adjust if necessary.
    for J = 1:ceil(m/B)
        I = 1+(J-1)*B : min(m,J*B);
        normAi(I) = sum(A(:,I).*A(:,I),1);
    end
elseif isa(A,'function_handle')
    for i = 1:m
        e = zeros(m,1); e(i) = 1;
        v = A(e,'transp');
        normAi(i) = norm(v)^2;
    end
else
    for i = 1:m
        e = zeros(m,1); e(i) = 1;
        v = A'*e;
        normAi(i) = norm(v)^2;
    end
end

% Depending on the ART method, set the row order.
is_randkaczmarz = false;
if ischar(art_method)
    switch lower(art_method)
        case 'kaczmarz'
            I = find(normAi>0);
        case 'symkaczmarz'
            I = find(normAi>0);
            I = [I, I(end:-1:1)];
        case 'randkaczmarz'
            is_randkaczmarz = true;
            I = find(normAi>0);
        otherwise
            error('Unknown ART method specified')
    end
else
    % Custom ART method specified by the user giving the row order I as
    % first input instead of the string with a particular ART method name.
    I = art_method(:)';
    I = I(normAi(I)>0);
end

% Apply damping.
normAi = normAi + damp*max(normAi);

% Special for randkaczmarz - set up random row selection.
if is_randkaczmarz
    cumul = cumsum(normAi/sum(normAi));
    
    % If all rows have approximately the same norm, treat them as having the
    % same norm when computing the random row index (which is much faster).
    if norm(cumul-(1:m)/m,inf) < 0.05  % Arbitrary threshold.
        fast = true;
    else
        fast = false;
    end
end

% Initalize waitbar if selected.
if do_waitbar
    h_waitbar = waitbar(0);
end

% Main ART loop.
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
    
    % Special for randkaczmarz - need random permutation or rows.
    if is_randkaczmarz
        I_torun = I(randperm(length(I)));
    else
        I_torun = I;
    end
    
    % The Kaczmarz sweep.
    for i = I_torun
        
        % Special for randkaczmarz.
        if is_randkaczmarz
            % The random row index.
            if fast
                ri = i;
            else
                ri = sum(cumul<rand)+1;
            end
        else
            % All other methods, just use i.
            ri = i;
        end
        
        if isnumeric(A)
            ai = A(:,ri); % Remember that A is transposed.
        elseif isa(A,'function_handle')
            e = zeros(m,1); e(ri) = 1;
            ai = Afun(e,'transp');  % ai is a column vector.
        else
            e = zeros(m,1); e(ri) = 1;
            ai = A'*e;  % ai is a column vector.
        end
        if isa(relaxparinput,'function_handle')
            relaxpar = relaxparinput((k-1)*m+i);
        end
        xk = xk + (relaxpar*(b(ri) - ai'*xk)/normAi(ri))*ai;
        
        % Enforce any lower and upper bounds (scalars or xk-sized vectors).
        if ~isempty(lbound)
            xk = max(xk,lbound);
        end
        if ~isnan(ubound)
            xk = min(xk,ubound);
        end
    end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
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

% Special for symkaczmarz: double finaliter.
if ischar(art_method) && strncmpi(art_method,'sym',3)
    info.finaliter = 2*info.finaliter;
    K = 2*K;
end

% List of iterates saved: all in K smaller than the final, and the final.
K = K(:); info.itersaved = [K(K<info.finaliter); info.finaliter];

% Save time total time taken.
info.timetaken = toc(t_total);