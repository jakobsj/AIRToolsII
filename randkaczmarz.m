function [X,info] = randkaczmarzPCH(A,b,K,x0,options)
%RANDKACZMARZ Randomized Kaczmarz method
%
%   [X,info] = randkaczmarz(A,b,K)
%   [X,info] = randkaczmarz(A,b,K,x0)
%   [X,info] = randkaczmarz(A,b,K,x0,options)
%
% Implements the randomized Kaczmarz's method for the system Ax = b:
%
%       x^{k+1} = x^k + lambda*(b_i - a^i*x^k)/(||a^i||_2^2)a^i
%
% where row i is chosen randomly with probability proportional to ||a^i||_2^2.
% One iteration consists of m such steps.
%
% Input:
%   A        m times n matrix, or a function that implements matrix-vector
%            multiplication with A and A'; please see explanation below.
%   b        m times 1 vector.
%   K        Number of iterations. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are saved, together with the last iterate.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       lambda    The relaxation parameter. For this method lambda must
%                 be a scalar; the default value is 1.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Peridogram.
%                            'DP'   : Discrepancy Principle.
%                     taudelta = the product of tau and delta, only
%                                necessary for DP.
%       nonneg    Logical; if true then nonnegativity in enforced in
%                 each step.
%       ubound    Upper bound in box constraint [0,ubound] on pixel values.
%       damping   A parameter D to avoid division by very small row norms
%                 by adding D*max_i{||a^i||_2^2} to ||a^i||_2^2.
%
% Output:
%   X     Matrix containing the saved iterations.
%   info  Information vector with 2 elements.
%         info(1) = 0 : stopped by maximum number of iterations
%                   1 : stopped by NCP-rule
%                   2 : stopped by DP-rule
%         info(2) = no. of iterations.
%         info(3) = the chosen lambda.
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
% 3) Then randkaczmarz is called with this A.
%
% See also: kaczmarz, symkaczmarz.

% Maria Saxild-Hansen and Per Chr. Hansen, Nov. 8, 2015, DTU Compute.

% Reference: T. Strohmer and R. Vershynin, A randomized Kaczmarz algorithm
% for linear systems with exponential convergence, J. Fourier Analysis and
% Applications, 15 (2009), pp. 262-278.

% Initialization.
if isa(A,'function_handle')
    mn = A([],'size');
    m = mn(1); n = mn(2);
else
    [m,n] = size(A);
    A = A';  % Faster to perform sparse column operations.
end

if nargin < 3
    error('Too few input arguments')
end

% Check that K is specified and initialize X matrix..
if isempty(K)
    error('Max no. of iterations must ALWAYS be specified')
end
Knew = sort(K);
kmax = Knew(end);
X = zeros(n,length(K));

% Default value of starting vector x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% The sizes of A, b and x must match.
if size(b,1) ~= m || size(b,2) ~= 1
    error('The sizes of A and b do not match')
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

% Initialization.
if nargin < 5

    stoprule = 'NO';
    lambda = 1;
    
    % Default is no nonnegativity or box constraint or damping.
    nonneg = false;
    boxcon = false;
    damp = 0;

else
% Check the contents of options, if present.
    
    % Nonnegativity.
    if isfield(options,'nonneg')
        nonneg = options.nonneg;
    else
        nonneg = false;
    end
    
    % Box constraints [0,L].
    if isfield(options,'ubound')
        nonneg = true;
        boxcon = true;
        L = options.ubound;
    else
        boxcon = false;
    end
    
    % Damping.
    if isfield(options,'damping')
        damp = options.damping;
        if damp<0, error('Damping must be positive'), end
    else
        damp = 0;
    end
    
    if isfield(options,'lambda')
        if isnumeric(options.lambda)
            lambda = options.lambda;
            if lambda <= 0 || lambda >= 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The lambda value is outside the interval (0,2)');
            end
        else
            error('lambda must be numeric')
        end
    else
        lambda = 1;
    end

    % Stopping rules
    if isfield(options,'stoprule') && isfield(options.stoprule,'type')
        stoprule = options.stoprule.type;
        if ischar(stoprule)
            if strncmpi(stoprule,'DP',2)
                % DP stopping rule
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error('The factor taudelta must be specified when using DP')
                end
                
                % Check that the first iteration should be performed:
                % rk = b - A'*x0;  % Remember that A is transposed.
                if isa(A,'function_handle')
                    r = b - A(x0,'notransp');
                else
                    r = b - A'*x0;  % Remember that A is transposed.
                end
                nrk = norm(r);
                if nrk <= taudelta
                    info = [2 0 NaN];
                    X = x0;
                    return
                end % end the DP-rule.
                
            elseif strncmpi(stoprule,'NC',2)
                % NCP stopping rule.
                dk = inf;
                q = floor(m/2);
                c_white = (1:q)'./q;
                K = [K max(K)+1];
                
            elseif strncmpi(stoprule,'NO',2)
                % No stopping rule.
                
            else
                error('The chosen stopping rule is not valid')
            end % end different stopping rules.
            
        else
            error('The stoprule type must be a string')
        end % end stoprule is a string.
        
    else
        % No stopping rule specified.
        stoprule = 'NO';
        
    end % end stoprule type specified.
    
end % end if nargin includes options.

% Calculate the norm of each row in A. This calculation can require a
% lot of memory. The commented lines can be used instead; they are
% slower, but use less memory!
if ~isa(A,'function_handle')
    normAi = full(abs(sum(A.*A,1)));  % Remember that A is transposed.
%   normAi = zeros(m,1);
%   for i = 1:m
%       ai = full(A(:,i));
%       normAi(i) = norm(ai)^2;
%   end
else
    normAi = zeros(1,m);
    for i = 1:m
        e = zeros(m,1);
        e(i) = 1;
        v = A(e,'transp');
        normAi(i) = norm(v)^2;
    end
end
I = find(normAi>0);
cumul = cumsum(normAi/sum(normAi));
normAi = normAi + damp*max(normAi);

% If all rows have approximately the same norm, treat them as having the
% same norm when computing the random row index (which is much faster).
if norm(cumul-(1:m)/m,inf) < 0.05  % Arbitrary threshold.
    fast = true;
else
    fast = false;
end

% Initialization before iterations.
xk = x0;
stop = 0;
k = 0;
l = 0;
klast = 0;

while ~stop
    k = k + 1;
    if strncmpi(stoprule,'NC',2)
        xkm1 = xk;
    end
    
    % The randomized Kaczmarz sweep; remember that A is transposed.
    for i = I(randperm(length(I)))
        % The random row index.
        if fast
            ri = i;
        else
            ri = sum(cumul<rand)+1;
        end
        % The updating step.
        if normAi(ri) > 0
            % xk = xk + (lambda*(b(ri) - A(:,ri)'*xk)/normAi(ri))*A(:,ri);
            if isa(A,'function_handle')
                e = zeros(m,1); e(ri) = 1;
                ai = A(e,'transp');  % ai is a column vector.
            else
                ai = A(:,ri); % Remember that A is transposed.
            end
            xk = xk + (lambda*(b(ri) - ai'*xk)/normAi(ri))*ai;
            if nonneg, xk = max(xk,0); end
            if boxcon, xk = min(xk,L); end
        end
    end
    
    % Stopping rules:
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        % nrk = norm(b - A'*xk);  % Remember that A is transposed.
        if isa(A,'function_handle')
            r = b - A(xk,'notransp');
        else
            r = b - A'*xk;  % Remember that A is transposed.
        end
        nrk = norm(r);
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [2 k lambda];
            else
                info = [0 k lambda];
            end
        end
        
    elseif strncmpi(stoprule,'NC',2)
        % NCP stopping rule.
        % rkh = fft(b - A'*xk);  % Remember that A is transposed.
        if isa(A,'function_handle')
            rkh = fft( b - A(xk,'notransp') );
        else
            rkh = fft(b - A'*xk);  % Remember that A is transposed.
        end
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end
        
        if dk < norm(c-c_white)
            stop = 1;
            if k ~= kmax
                info = [1 k-1 lambda];
            else
                info = [0 k-1 lambda];
            end
        else
            dk = norm(c-c_white);
        end % end NCP-rule.
        
    elseif strncmpi(stoprule,'NO',2)
        % No stopping rule.
        if k >= kmax 
            stop = 1;
            info = [0 k lambda];
        end
    end % end stoprule type.
    
    % If the current iteration is requested saved.
    if k == Knew(l+1) || stop
        l = l + 1;
        % Saves the current iteration.
        if strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xkm1;
            else
                l = l - 1;
            end
        else
            X(:,l) = xk;
        end
        klast = k;
    end
    
end
X = X(:,1:l);