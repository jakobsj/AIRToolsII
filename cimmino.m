function [X,info,restart] = cimmino(A,b,K,x0,options)
%CIMMINO Cimmino's method
%
%   [X,info,restart] = cimmino(A,b,K)
%   [X,info,restart] = cimmino(A,b,K,x0)
%   [X,info,restart] = cimmino(A,b,K,x0,options)
%
% Implements Cimmino's method for the linear system Ax = b:
%
%       x^{k+1} = x^k + lambda_k*A'*M*(b-A*x^k) ,
%
% where M = (1/m)*diag(w_i/||a^i||_2^2) and w_i are weights (default w_i = 1).
%
% Input:
%   A        m times n matrix, or a function that implements matrix-vector
%            multiplication with A and A'; please see explanation below.
%   b        m times 1 vector containing the right-hand side.
%   K        Number of iterations. If K is a scalar, then K is the maximum
%            number of iterations and only the last iterate is saved.
%            If K is a vector, then the largest value in K is the maximum
%            number of iterations and all iterates corresponding to the
%            values in K are saved, together with the last iterate.no
%   x0       n times 1 starting vector. Default: x0 = 0.
%   options  Struct with the following fields:
%       lambda    The relaxation parameter. If lambda is a scalar then
%                 the corresponding value is used in each iteration;
%                 default value is 1.9/norm(A'*M*A).
%                 If lambda is a string, then it refers to a method to 
%                 determine lambda in each iteration. For this method the
%                 following strings can be specified:
%                     'line'    : lambda is chosen using line search.
%                     'psi1'    : lambda is chosen using the Psi_1-based 
%                                 relaxation method.
%                     'psi1mod' : lambda is chosen using the modified     
%                                 Psi_1-based relaxation method.
%                     'psi2'    : lambda is chosen using the Psi_2-based
%                                 relaxation method.
%                     'psi2mod' : lambda is chosen using the modified     
%                                 Psi_2-based relaxation method.
%       stoprule  Struct containing the following information about the
%                 stopping rule:
%                     type = 'none' : (Default) the only stopping rule
%                                     is the maximum number of iterations.
%                            'NCP'  : Normalized Cumulative Perodogram.
%                            'DP'   : Discrepancy Principle.
%                            'ME'   : Monotone Error rule.
%                     taudelta = product of tau and delta, only needed
%                                for DP and ME.
%       nonneg    Logical; if true then nonnegativity in enforced in
%                 each iteration.
%       ubound    Upper bound in box constraint [0,ubound] on pixel values.
%       restart   Struct that contains the first singular value s1 of
%                 sqrt(M)*A and the diagonals of the matrix M.
%       w         m-dimensional weigthing vector.
%
% Output:
%   X        Matrix containing the saved iterations.
%   info     Information vector with 2 elements.
%            info(1) = 0 : stopped by maximum number of iterations
%                      1 : stopped by NCP-rule
%                      2 : stopped by DP-rule
%                      3 : stopped by ME-rule.
%            info(2) = no. of iterations.
%            info(3) = the chosen lambda.
%   restart  Struct containing the largest singular value s1 and the
%            diagonal of the matrix M.
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
% 3) Then cimmino is called with this A.
%
% See also: landweber, cav, drop, sart.

% Maria Saxild-Hansen, Per Chr. Hansen and Jakob Sauer J�rgensen,
% November 8, 2015, DTU Compute.

% Reference: C. D. Meyer, Matrix Analysis and Applied Linear Algebra, SIAM,
% Philadelphia, 2000 (this reference uses reflections - instead of projec-
% tions and therefore includes a redundant factor of 2).

% Check that at least 3 inputs are given
if nargin < 3
    error('Too few input arguments')
end

% If A is not a function (i.e., A is a matrix or an object), convert A
% to a function.
if isa(A,'function_handle')
    Afun = A;
else
    Afun = @(xx,transp_flag) Afun_matrix(xx,transp_flag,A);
end

% Check that the sizes of A and b match.
mn = Afun([],'size');
m = mn(1); n = mn(2);
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
end

% Check that K is specified and initialize X matrix..
if isempty(K)
    error('Max no. of iterations must ALWAYS be specified')
end
Knew = sort(K);
kmax = Knew(end);
X = zeros(n,length(K));

% Default value for x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% Check the size of x0.
if size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem.')
end

rxk = b - Afun(x0,'notransp');

if nargin < 5
    
    % Calculate the norm of each row in A. This calculation can require a
    % lot of memory. The commented lines can be used instead; they are
    % slower, but use less memory!
    if ~isa(A,'function_handle')
        normAi = full(abs(sum(A.*A,2)));
%       normAi = zeros(m,1);
%       for i = 1:m
%           ai = full(A(i,:));
%           normAi(i) = norm(ai)^2;
%       end
    else
        normAi = zeros(m,1);
        for i = 1:m
            e = zeros(m,1);
            e(i) = 1;
            v = Afun(e,'transp');
            normAi(i) = norm(v)^2;
        end
    end
    
    % Define the M matrix.
    M = (1/m)*(1./normAi);
    I = (M == Inf);
    M(I) = 0;
    
    % Calculate the largest singular value.
    atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
    optionsEIGS.disp = 0;
    sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
    
    % If restart is required as output.
    if nargout == 3
        restart.M = M;
        restart.s1 = sqrt(sigma1tildesquare);
    end
    
    % Default value of lambda.
    lambda = 1.9/sigma1tildesquare;
    casel = 1;
    
    % Default stopping rule.
    stoprule = 'NO';
    k = 0;
    
    % Default is no nonnegativity or box constraint.
    nonneg = false;
    boxcon = false;

else
% Check the contents of options if present.
    
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
    
    if isfield(options,'restart') && isfield(options.restart,'M')
        M = options.restart.M;
    else
        % Calculate the norm of each row in A. This calculation can require 
        % a lot of memory. The commented lines can be used instead; they 
        % are slower, but use less memory!
        if ~isa(A,'function_handle')
            normAi = full(abs(sum(A.*A,2)));
%           normAi = zeros(m,1);
%           for i = 1:m
%               ai = full(A(i,:));
%               normAi(i) = norm(ai)^2;
%           end
        else
            normAi = zeros(m,1);
            for i = 1:m
                e = zeros(m,1);
                e(i) = 1;
                v = Afun(e,'transp');
                normAi(i) = norm(v)^2;
            end
        end
        
        % If the method is weigthed.
        if isfield(options,'w')
            w = options.w;
            M = 1/m*(w./normAi);
        else
            M = 1/m*(1./normAi);
        end
        
        % Replace the Inf in the zero rows with 0.
        I = (M == Inf);
        M(I) = 0;
    end
    
    % If restart is required as output.
    if nargout == 3
        restart.M = M;
    end    
    
    % Stopping rules
    if isfield(options,'stoprule') && isfield(options.stoprule,'type')
        stoprule = options.stoprule.type;
        % Check that the stoprule is a string.
        if ischar(stoprule)
            if strncmpi(stoprule,'DP',2)
                % DP stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error(['The factor taudelta must be specified ',...
                        'when using DP'])
                end
                k = 0;
                
                % Check that the first iteration should be performed:
                rk = rxk;
                nrk = norm(rk);
                
                if nrk <= taudelta
                    info = [2 k NaN];
                    X = x0;
                    return
                end % end the DP-rule.
                
            elseif strncmpi(stoprule,'ME',2)
                % ME stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error(['The factor taudelta must be specified ',...
                        'when using ME'])
                end
                k = 0;
                rk = rxk;
                K = K + 1;
                
            elseif strncmpi(stoprule,'NC',2)
                % NCP stopping rule.
                dk = inf;
                q = floor(m/2);
                c_white = (1:q)'./q;
                k = 0;
                K = [K  max(K)+1];
                
            elseif strncmpi(stoprule,'NO',2)
                k = 0;
                
            else
                error('The chosen stopping rule is not valid.')
            end % end stoprule is string.
            
        else
            error('The stoprule type must be a string')
        end % end stoprule is string.
        
    else
        % No stopping rule specified.
        stoprule = 'NO';
        k = 0;
              
    end % end stoprule type specified.
          
    % Determine the relaxation parameter lambda.
    if isfield(options,'lambda')
        lambda = options.lambda;
        
        % If lambda is a scalar.
        if ~ischar(lambda)
            % Check if the first singular value is known in options.
            if isfield(options,'restart') && isfield(options.restart,'s1')
                sigma1tildesquare = options.restart.s1^2;
            else
                % Calculate the first singular value.
                atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
                optionsEIGS.disp = 0;
                sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
            end
            
            % If restart is required as output.
            if nargout == 3
                restart.s1 = sqrt(sigma1tildesquare);
            end
            
            % Convergence check.
            if lambda <= 0 || lambda >= 2/sigma1tildesquare;
                warning('MATLAB:UnstableRelaxParam',...
                    ['The lambda value is outside the ',...
                    'interval (0,%f)'],2/sigma1tildesquare);
            end
            casel = 1;
        else
            % Calculates the lambda value according to the chosen method.
            if strncmpi(lambda,'line',4)
                % Method: Line search.
                casel = 2;
                
                if nargout == 3
                    if isfield(options,'restart') && ...
                            isfield(options.restart,'s1')
                        sigma1tildesquare = options.restart.s1^2;
                    else
                        % Calculate the first singular value.
                        atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
                        optionsEIGS.disp = 0;
                        sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
                    end
                    restart.s1 = sqrt(sigma1tildesquare);
                end
            elseif strncmpi(lambda,'psi1',4)
                % Method: ENH psi1.
                casel = 3;
                ite = 0;
                
                % Checks if the largest singular value is known.
                if isfield(options,'restart') && ...
                        isfield(options.restart,'s1')
                    sigma1tildesquare = options.restart.s1^2;
                else
                    % Calculate the first singular value.
                    atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
                    optionsEIGS.disp = 0;
                    sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
                end
                
                % If restart is required as output.
                if nargout == 3
                    restart.s1 = sqrt(sigma1tildesquare);
                end
                
                % Precalculate the roots.
                z = calczeta(2:max(K)-1);
                
                % Define the values for lambda according to the psi1
                % strategy modified or not.
                if strncmpi(lambda,'psi1mod',7)
                    nu = 2;
                    lambdak = [sqrt(2); sqrt(2); 
                                nu*2*(1-z)]/sigma1tildesquare;
                else
                    lambdak = [sqrt(2); sqrt(2); 
                                2*(1-z)]/sigma1tildesquare;
                end
                
            elseif strncmpi(lambda,'psi2',4)
                % Method: ENH psi2.
                casel = 3;
                ite = 0;
                
                % Checks if the largest singular value is known.
                if isfield(options,'restart') && ...
                        isfield(options.restart,'s1')
                    sigma1tildesquare = options.restart.s1^2;
                else
                    % Calculate the first singular value.
                    atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
                    optionsEIGS.disp = 0;
                    sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
                end
                
                % If restart is required as output.
                if nargout == 3
                    restart.s1 = sqrt(sigma1tildesquare);
                end
                
                % Precalculate the roots.
                kk = 2:max(K)-1;
                z = calczeta(kk);
                
                % Define the values for lambda according to the psi2
                % strategy modified or not.
                if strncmpi(lambda,'psi2Mod',7)
                    nu = 1.5;
                    lambdak = [sqrt(2); sqrt(2); 
                           nu*2*(1-z)./((1-z.^(kk')).^2)]/...
                           sigma1tildesquare;
                else
                    lambdak = [sqrt(2); sqrt(2);
                           2*(1-z)./((1-z.^(kk')).^2)]/sigma1tildesquare;
                end
                
            else
                error(['The chosen relaxation strategy is not valid '...
                    'for this method.'])                      
            end % end check of lambda strategies.
            
        end % end check of the class of lambda.
    else
        % Check if the first singular value is given.
        if isfield(options,'restart') && isfield(options.restart,'s1')
            sigma1tildesquare = options.restart.s1^2;
        else
            % Calculates the first singular value
            atma = @(x) Afun( M.*Afun(x,'notransp') , 'transp' );
            optionsEIGS.disp = 0;
            sigma1tildesquare = eigs(atma,n,1,'lm',optionsEIGS);
        end
        
        % If restart is required as output.
        if nargout == 3
            restart.s1 = sqrt(sigma1tildesquare);
        end
        
        % Define a default constant lambda value.
        lambda = 1.9/sigma1tildesquare;
        casel = 1;
        
    end % end if lambda is a field in options.

end % end if nargin includes options.

% Initialize the values
xk = x0;
stop = 0;
l = 0;
klast = 0;

while ~stop
    % Update the iteration number k.
    k = k + 1;
    if strncmpi(stoprule,'ME',2)  || strncmpi(stoprule,'NC',2)
        xkm1 = xk;
    end
    % Compute the current iteration.
    if casel == 1
        % Cimmino using constant value of lambda.
        xk = xk + lambda*Afun(M.*rxk,'transp'); %(A'*(M.*rxk));
    elseif casel == 2
        % Cimmino using line search.
        Mrk = M.*rxk;
        ATMrk = Afun(Mrk,'transp'); %A'*Mrk;
        lambdak = (rxk'*Mrk)/norm(ATMrk)^2;
        xk = xk + lambdak*ATMrk;
    elseif casel == 3
        % Cimino using psi1 og psi2.
        ite = ite + 1;
        xk = xk + lambdak(ite)*Afun(M.*rsk,'transp'); %(A'*(M.*rxk));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rxk = b - Afun(xk,'notransp'); %A*xk;
    
    % Stopping rules:
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        nrk = norm(rxk);
                
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [2 k lambda];
            else
                info = [0 k lambda];
            end
        end % end the DP-rule.
        
    elseif strncmpi(stoprule,'ME',2)
        % ME stopping rule.
        nrk = norm(rk);
        dME = rk'*(rk+rxk)/2;
    
        if dME/nrk <= taudelta || k >= kmax
            stop = 1;
            
            if k ~= kmax
                info = [3 k-1 lambda];
            else
                info = [0 k-1 lambda];
            end        
        else
            rk = rxk;
        end % end the ME-rule.
        
    elseif strncmpi(stoprule,'NC',2)
        % NCP stopping rule.
        rkh = fft(rxk);
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end

        if dk < (norm(c-c_white)) || k >= kmax
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
        % Save the current iteration.
        if strncmpi(stoprule,'ME',2)
                X(:,l) = xkm1;

        elseif strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xkm1;
            else
                % Since the iteration was not saved.
                l = l - 1;
            end
        else
            X(:,l) = xk;
        end
        klast = k;
    end
    
end
X = X(:,1:l);