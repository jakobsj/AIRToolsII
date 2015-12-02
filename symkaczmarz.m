function [X,info] = symkaczmarzPCH(A,b,K,x0,options)
%SYMKACZMARZ Symmetric Kaczmarz method
%
%   [X,info] = symkaczmarz(A,b,K)
%   [X,info] = symkaczmarz(A,b,K,x0)
%   [X,info] = symkaczmarz(A,b,K,x0,options)
%
% Implements the symmetric Kaczmarz method: each sweep consists of a Kaczmarz
% sweep followed by a Kaczmarz sweep with the equations in reverse order.
%
% Input:
%   A        m times n matrix, or a function that implements matrix-vector
%            multiplication with A and A'; please see explanation below.
%   b        m times 1 vector.
%   K        Number of iteration. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and only iterates corresponding to the
%            values in K are saved, together with the last iterate.
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       lambda    The relaxation parameter. If lambda is a scalar then
%                 the corresponding value is used in each iteration; the
%                 default value is 1.
%                 If lambda is a string, then it refers to a method to 
%                 determine lambda in each iteration. For this method the
%                 following strings can be specified:
%                     'psi1' : lambda is chosen using the Psi_1-based 
%                                 relaxation method.
%                     'psi2' : lambda is chosen using the Psi_2-based
%                                 relaxation method.
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
% 3) Then symkaczmarz is called with this A.
%
% See also: kaczmarz, randkaczmarz

% Maria Saxild-Hansen and Per Chr. Hansen, Nov. 8, 2015, DTU Compute.

% Reference: Å. Björck and T. Elfving, Accelerated projection methods for 
% computiong pseudoinverse solutions of systems of linear equations, BIT 
% 19 (1979), pp. 145-163.

% Initialization.
if isa(A,'function_handle')
    mn = A([],'size');
    m = mn(1); n = mn(2);
else
    [m,n] = size(A);
    A = A';  % Faster to perform sparse column operations.
end

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

% Input check: The sizes of A, b and x must match.
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

if nargin < 5

    stoprule = 'NO';
    lambda = 1;
    casel = 1;
    
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
   
    % Stopping rules.
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
                error('The shosen stopping rule is not valid')
            end % end different stopping rules
            
        else
            error('The stoprule type must be a string')
        end % end stoprule is a string.

    else
        % No stopping rule specified.
        stoprule = 'NO';
        
    end % end stoprule type specified.
    
    if isfield(options,'lambda')
        lambda = options.lambda;
        % If lambda is a scalar.
        if ~ischar(lambda)
            % Convergence check.
            if lambda < 0 || lambda > 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The lambda value is outside the interval (0,2)');
            end
            casel = 1;
        else
            % Calculates the lambda value according to the chosen method.
            if strncmpi(lambda,'psi1',4)
                % Method: Psi1
                casel = 3;
                ite = 0;
                
                sigma1tildesquare = 1;
                       
                % Precalculates the roots.
                z = calczeta(2:max(K)-1);
                
                % Define the values for lambda according to Psi1.
                lambdak = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tildesquare;
                
            elseif strncmpi(lambda,'psi2',4)
                % Method: Psi2.
                casel = 3;
                ite = 0;
                
                sigma1tildesquare = 1;
                
                % Precalculates the roots.
                kk = 2:max(K)-1;
                z = calczeta(kk);
                
                % Define the values for lambda according to the psi2.
                lambdak = [sqrt(2); sqrt(2); 
                    2*(1-z)./((1-z.^(kk')).^2)]/sigma1tildesquare;
            else
                error(['The chosen relaxation strategy is not ',...
                    'valid for this metod.'])
            end % end check of lambda strategies.
        end
        
    else
        casel = 1;
        lambda = 1;
    end
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
I = [I, I(end-1:-1:2)];
normAi = normAi + damp*max(normAi);

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
    
    % The Kaczmarz sweep followed by the reverse order sweep.
    for i = I
        if isa(A,'function_handle')
            e = zeros(m,1); e(i) = 1;
            ai = A(e,'transp');  % ai is a column vector.
        else
            ai = A(:,i); % Remember that A is transposed.
        end
        if casel == 1
            xk = xk + (lambda*(b(i) - ai'*xk)/normAi(i))*ai;
        else
            ite = ite + 1;
            xk = xk + (lambdak(ite)*(b(i) - ai'*xk)/normAi(i))*ai;
        end
        if nonneg, xk = max(xk,0); end
        if boxcon, xk = min(xk,L); end
    end
    
    % Stopping rules:
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        % nrk = norm(b-A'*xk);  % Remember that A is transposed.
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
        % rkh = fft( b-A'*xk);  % Remember that A is transposed.
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
        
        if dk < norm(c-c_white) || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [1 k-1 lambda];
            else
                info = [0 k-1 lambda];
            end
        else
            dk = norm(c-c_white);
        end
        
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
        % Savs the current iteration.
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