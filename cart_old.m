function [X,info] = cart_old(A,b,K,x0,options)
%CART Columnwise version of Kaczmarz's method
%
%   [X,info] = cart(A,b,K)
%   [X,info] = cart(A,b,K,x0)
%   [X,info] = cart(A,b,K,x0,options)
%
% Implements the CART method for the system Ax = b:
%
%       x_j^{k+1} = x_j^k + lambda*a^j'(b - A x^k)/(||a^j||_2^2)
%
% where a^j is the j-th column of A, and j = (k mod n) + 1.
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
%       lambda    The relaxation parameter. For this method lambda must
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
%         info(3) = the chosen lambda.
%
% See also: kaczmarz, randkaczmarz, symkaczmarz

% Jacob Fr�sig, Nicolai Riis, Per Chr. Hansen, Nov. 7, 2015, DTU Compute.

if isa(A,'function_handle'), error('A cannot be a function handle'), end
[m,n] = size(A);

% Check the number of inputs.
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
    error('The size of A and b do not match')
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

% Initialization.
if nargin < 5

    stoprule = 'NO';
    lambda = 0.25;
    
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
        lambda = 0.25;
    end
    
    % Stopping rules
    if isfield(options,'stoprule') && isfield(options.stoprule,'type')
        stoprule = options.stoprule.type;
        if ischar(stoprule)
            if strncmpi(stoprule,'DP',2)
                % DP stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error('The factor taudelta must be specified when using DP')
                end
                        
                % Check that the first iteration should be performed:
                rk = (b-A*x0);
                nrk = norm(rk);
                
                if nrk <= taudelta
                    info = [2 0 NaN];
                    X = x0;
                    return
                end % end the DP-rule.
                
            elseif strncmpi(stoprule,'NC',2)
                % NCP stopping rule.
                dk = inf;
                q = floor(m/2);
                c_white = (1:q)./q;
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

% Initialization before iterations.
xk = x0;
normAj = full(abs(sum(A.*A,1)));
normAj = normAj + damp*max(normAj);
J = find(normAj>0);

stop = 0;
k = 0;
l = 0;
klast = 0;
r = b - A*x0;

while ~stop
    k = k + 1;
    if strncmpi(stoprule,'NC',2)
        xkm1 = xk;
    end
    % The columnwise Kaczmarz sweep.
    for j = J
        delta = A(:,j)'*r/normAj(j);
        od = lambda*delta;
        
        if nonneg && od < -xk(j)
            od = -xk(j);
        end
        if boxcon && od > L - xk(j)
            od = L - xk(j);
        end
        xk(j) = xk(j) + od;
        r = r - od*A(:,j);
    end
    
    % Stopping rules.
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        nrk = norm(b - A*xk);
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k >= kmax
                info = [0 k lambda];
            else
                info = [2 k lambda];
            end
        end % end the DP-rule.
        
    elseif strncmpi(stoprule,'NC',2)
        % NCP stopping rule.
        rkh = fft(b - A*xk);
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end
        if dk < norm(c-c_white') || k >= kmax
            stop = 1;
            if k >= kmax
                info = [0 k-1 lambda];
            else
                info = [1 k-1 lambda];
            end
        else
            dk = norm(c-c_white');
        end % end NCP-rule.
        
    elseif strncmpi(stoprule,'NO',2)
        % No stopping rule.
        if k >= kmax
            stop = 1;
            info = [0 k lambda];
        end
    end % end stoprule type.
    
    % If the current iteration is requested saved.
    if (~isempty(K) && k == Knew(l+1)) || stop
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