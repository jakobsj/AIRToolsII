function [X,info,restart] = landweber(A,b,K,x0,options)
%LANDWEBER The classical Landweber method
%
%   [X,info,restart] = landweber(A,b,K)
%   [X,info,restart] = landweber(A,b,K,x0)
%   [X,info,restart] = landweber(A,b,K,x0,options)
%
% Implements the classical Landweber method for the linear system Ax = b:
%
%       x^{k+1} = x^k + lambda_k*A'*(b-A*x^k)
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
%       lambda    The relaxation parameter. If lambda is a scalar then
%                 the corresponding value is used in each iteration;
%                 default value is 1.9/norm(A)^2.
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
%       restart   Struct that can contain the largest singular value s1 of A.
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
%   restart  Struct containing the largest singular value s1.
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
% 3) Then landweber is called with this A.
%
% See also: cimmino, cav, drop, sart.

% Maria Saxild-Hansen, Per Chr. Hansen and Jakob Sauer J�rgensen,
% November 8, 2015, DTU Compute.

% Reference: L. Landweber, An iteration formula for Fredholm integral
% equations of the first kind, American Journal of Mathematics, 73 (1951),
% pp. 615-624.

% Check inputs.
[Afun,m,n,Knew,kmax,x0] = check_inputs(A,b,K,x0,options);

X = zeros(n,length(Knew));

rxk = b - Afun(x0,'notransp');

if nargin < 5
    
    % Calculate the largest singular value.
    atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
    optionsEIGS.disp = 0;
    sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
    
    % If restart is required as output.
    if nargout == 3
        restart.s1 = sigma1tilde;
    end
    
    % Default value of lambda.
    lambda = 1.9/sigma1tilde^2;
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
    
    % Initialize for stopping rules.
    [k,K,rk,dk] = init_stoprules(stoprule,K);
    
    % Do initial check of stopping criteria - probably lambda should be set
    % before this, perhaps just to nan.
    [stop, info, rk, dk] = check_stoprules(...
    stoprule, rxk, lambda, taudelta, k, kmax, rk, dk)
    
    % Determine the relaxation parameter lambda.
    if isfield(options,'lambda')
        lambda = options.lambda;
        
        % If lambda is a scalar.
        if ~ischar(lambda)
            % Checks if the largest singular value is known in options.
            if isfield(options,'restart') && isfield(options.restart,'s1')
                sigma1tilde = options.restart.s1;
            else
                % Calculates the largest singular value.
                % sigma1tilde = normest_mf(Afun);
                atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
                optionsEIGS.disp = 0;
                sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
            end
            
            % If restart is required as output.
            if nargout == 3
                restart.s1 = sigma1tilde;
            end
            
            % Checks if the given constant lambde value is unstable.
            if lambda <= 0 || lambda >= 2/sigma1tilde^2
                warning('MATLAB:UnstableRelaxParam',['The lambda value '...
                    'is outside the interval (0,%f)'],2/sigma1tilde^2)
            end
            casel = 1;
        else
            % Calculate the lambda value according to the chosen method.
            if strncmpi(lambda,'line',4)
                % Method: Line search
                casel = 2;
                
                if nargout == 3
                    if isfield(options,'restart') && ...
                            isfield(options.restart,'s1')
                        sigma1tilde = options.restart.s1;
                    else
                        % sigma1tilde = normest_mf(Afun);
                        atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
                        optionsEIGS.disp = 0;
                        sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
                        restart.s1 = sigma1tilde;
                    end
                    restart.s1 = sigma1tilde;
                end
                
            elseif strncmpi(lambda,'psi1',4)
                % Method: ENH psi1.
                casel = 3;
                ite = 0;
                
                % Checks if the largest singular value is known.
                if isfield(options,'restart') && ...
                        isfield(options.restart,'s1')
                    sigma1tilde = options.restart.s1;
                else
                    % Calculates the largest singular value.
                    % sigma1tilde = normest_mf(Afun);
                    atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
                    optionsEIGS.disp = 0;
                    sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
                end
                
                % If restart is required as output.
                if nargout == 3
                    restart.s1 = sigma1tilde;
                end
                
                % Precalculate the roots.
                z = calczeta(2:max(K)-1);
                
                % Define the values for lambda according to the psi1
                % strategy modified or not.
                if strncmpi(lambda,'psi1mod',7)
                    nu = 2;
                    lambdak = [sqrt(2); sqrt(2); nu*2*(1-z)]/sigma1tilde^2;
                else
                    
                    lambdak = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tilde^2;
                end
                
            elseif strncmpi(lambda,'psi2',4)
                % Method: ENH psi2.
                casel = 3;
                ite = 0;
                
                % Checks if the largest singular value is known.
                if isfield(options,'restart') && ...
                        isfield(options.restart,'s1')
                    sigma1tilde = options.restart.s1;
                else
                    % Calculates the largest singular value.
                    % sigma1tilde = normest_mf(Afun);
                    atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
                    optionsEIGS.disp = 0;
                    sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
                end
                
                % If restart is required as output.
                if nargout == 3
                    restart.s1 = sigma1tilde;
                end
                
                % Precalculate the roots.
                kk = 2:max(K)-1;
                z = calczeta(kk);
                
                % Define the values for lambda according to the psi2
                % strategy modified or not.
                if strncmpi(lambda,'psi2Mod',7)
                    nu = 1.5;
                    lambdak = [sqrt(2); sqrt(2);
                        nu*2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
                else
                    lambdak = [sqrt(2); sqrt(2);
                        2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
                end
                
            else
                error(['The chosen relaxation strategy is not valid '...
                    'for this method.'])
            end % end check of the class of lambda.
            
        end % end check of lambda strategies.
    else
        % Check if the largest singular value is given.
        if isfield(options,'restart') && isfield(options.restart,'s1')
            sigma1tilde = options.restart.s1;
        else
            % Calculates the largest singular value.
            % sigma1tilde = normest_mf(Afun);
            atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
            optionsEIGS.disp = 0;
            sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
        end
        
        % If restart is required as output.
        if nargout == 3
            restart.s1 = sigma1tilde;
        end
        
        % Define a default constant lambda value.
        lambda = 1.9/sigma1tilde^2;
        casel = 1;
        
    end % end if lambda is a field in options.
    
end % end if nargin includes options.

% Initialize the values.
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
        % Landweber using constant value of lambda.
        xk = xk + lambda*(Afun(rxk,'transp'));
    elseif casel == 2
        % Landweber using line search.
        ATrk = Afun(rxk,'transp');
        lambdak = norm(rxk)^2/norm(ATrk)^2;
        xk = xk + lambdak*ATrk;
    elseif casel == 3
        % Landweber using psi1 or psi2.
        ite = ite + 1;
        xk = xk + lambdak(ite)*(Afun(rxk,'transp'));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rxk = b - Afun(xk,'notransp');
    
    % Check stopping rules:
    [stop,info,rk,dk] = check_stoprules(...
        stoprule, rxk, lambda, taudelta, k, kmax, rk, dk);
        
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