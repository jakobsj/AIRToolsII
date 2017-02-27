function [X,info] = cart(varargin)
%CART Columnwise version of Kaczmarz's method
%
%   [X,info] = cart(A,b,K)
%   [X,info] = cart(A,b,K,x0)
%   [X,info] = cart(A,b,K,x0,options)
%
% Implements the CART method for the system Ax = b:
%
%       x_j^{k+1} = x_j^k + relaxpar*a^j'(b - A x^k)/(||a^j||_2^2)
%
% where a^j is the j-th column of A, and j = (k mod n) + 1.
%
% This version incorporates the "flagging" idea: a component of the
% solution x is "flagged" if its update is smaller than a threshold THR
% times max(abs(x)) is not updated again until it is "unflagged" which
% happens after a random number of iterations chosen as rand*Nunflag.
%
% Input:
%   A        m times n matrix.
%   b        m times 1 vector.
%   K        Number of iterations. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are saved, together with the last iterate.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       relaxpar    The relaxation parameter. For this method relaxpar must
%                 be a scalar; default value is 0.25.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Peridogram.
%                            'DP'   : Discrepancy Principle.
%                       taudelta = The product of tau and delta, only
%                                  necessary for DP.
%       nonneg    Logical; if true then nonnegativity in enforced in
%                 each step.
%       ubound    Upper bound in box constraint [0,ubound] on pixel values.
%       damping   A parameter D to avoid division by very small column norms
%                 by adding D*max_j{||a^j||_2^2} to ||a^j||_2^2.
%
% Output:
%   X     Matrix containing the saved iterations.
%   info  Information vector with 2 elements.
%         info(1) = 0 : stopped by maximum number of iterations
%                   1 : stopped by NCP-rule
%                   2 : stopped by DP-rule
%         info(2) = no. of iterations.
%         info(3) = the chosen relaxpar.
%
% See also: kaczmarz, randkaczmarz, symkaczmarz

% Jacob Fr�sig, Nicolai Riis, Per Chr. Hansen, Nov. 7, 2015, DTU Compute.



% TODOS
% 
% - In cartFlag third output of info is NOT relaxpar, but "number of times an
% element is updated."
% 
% - In cartFlag additional options.THR=1e-4, options.Kbegin=10,
% options.Nunflag=max(K)/4
% 
% - In old kaczmarz, I is computed before applying damping. In cart J is
% only computed AFTER damping is applied. Not consistent? Correct?
% 
% - Implement vector lbound/ubound


% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    s1,w,res_dims,ncp_smooth,damp,THR,Kbegin,Nunflag] = ...
    check_inputs(varargin{:});

A = varargin{1};

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
if strcmpi(stoprule,'ME')
    error('Stopping rule ME is not available for CART method.')
end
[k,rkm1,dk] = init_stoprules(stoprule,rk,ncp_smooth);

% Do initial check of stopping criteria - probably relaxpar should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Compute the relaxation parameter to be used throughout iterations.
relaxpar = calcrelaxpar_cart(relaxparinput);

% Calculate the norm of each column in A. This calculation can require a
% lot of memory. Unlike art methods, A is NOT transposed.
if ~isa(A,'function_handle')
    normAj = full(abs(sum(A.*A,1)));
else
    normAj = zeros(1,n);
    for j = 1:n
        e = zeros(n,1);
        e(j) = 1;
        v = A(e,'notransp');
        normAj(j) = norm(v)^2;
    end
end

% Apply damping and determine non-zero columns.
normAj = normAj + damp*max(normAj);
J = find(normAj>0);

% Initialization before iterations.
xk = x0;
l = 1;

% Set additional CART "flagging" parameters.
F = true(n,1);       % Vector of logical "flags."
Nflag = zeros(n,1);  % Counts "flagged" iterations for each component.
DOT = 0;
AXPY = 0;
UPD = 0;
SKIP = 0;
DODO = zeros(n,max(K));
WORK = zeros(max(K),1);


while ~stop
    
    % Update the iteration number k.
    k = k + 1;
    
    % The columnwise Kaczmarz sweep.
    for j = J
        mm = max(abs(xk)); % Max abs element of previous iteration.
        if F(j)
            DOT = DOT + 1;
            
            % Get the j'th column of A
            if ~isa(A,'function_handle')
                aj = A(:,j);
            else
                e = zeros(n,1);
                e(j) = 1;
                aj = A(e,'notransp');
            end
            
            delta = aj'*rk/normAj(j);
            od = relaxpar*delta;           % The update.
        
            % Correction for constraints.
            if ~isnan(lbound) && od < lbound - xk(j)    
                od = lbound - xk(j);
            end
            if ~isnan(ubound) && od > ubound - xk(j)
                od = ubound - xk(j);
            end
            xk(j) = xk(j) + od;
            
            DODO(j,k) = abs(od);
            UPD = UPD + 1;
            
            if k >= Kbegin && abs(od) < THR*mm; % "Flag" if needed.
                F(j) = false;
                Nflag(j) = 1;  % PCH
            end

            
            AXPY = AXPY + 1;
            rk = rk - od*aj;
        else
            SKIP = SKIP + 1;
            if Nflag(j) < randi(Nunflag)
                Nflag(j) = Nflag(j) + 1;
            else
                F(j) = true;
            end
        end
    end
    WORK(k) = DOT + AXPY;
    
    % New residual.
    %rk = b - Afun(xk,'notransp');
    %rk = b - Afun(xk,'notransp');
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims);
    
    % If the current iteration is requested saved.
    if k == K(l) || stop
        X(:,l) = xk;
        l = l + 1;
    end
end

% Return only the saved iterations: Only to "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);

% Special for CART: Update info.
info.DOT = DOT;
info.AXPY = AXPY;
info.UPD = UPD;
info.SKIP = SKIP;