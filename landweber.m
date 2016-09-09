function [X,info,restart] = landweber(varargin)
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

% Maria Saxild-Hansen, Per Chr. Hansen and Jakob Sauer Joergensen,
% November 8, 2015, DTU Compute.

% Reference: L. Landweber, An iteration formula for Fredholm integral
% equations of the first kind, American Journal of Mathematics, 73 (1951),
% pp. 615-624.

% Parse inputs.
[Afun,b,m,n,K,Knew,kmax,x0,nonneg,boxcon,L,stoprule,taudelta,...
    lambdainput,s1,dims] = check_inputs(varargin{:});

X = zeros(n,length(Knew));

% TODO Should there be both Knew and K, or reduce to just Knew?

rxk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
[k,K,rk,dk] = init_stoprules(stoprule,rxk,K);

% Do initial check of stopping criteria - probably lambda should be set
% before this, perhaps just to nan.
[stop, info, rk, dk] = check_stoprules(...
    stoprule, rxk, lambdainput, taudelta, k, kmax, rk, dk, dims);

% TODO If is aborting here, is output X set? Perhaps make sure x0 is always
% written to the first column of X.

% Calculate lambda and restart
atma = @(x) Afun( Afun(x,'notransp') , 'transp' );
[lambda, casel, restart] = calclambda(lambdainput, s1, K, atma, n);

% Initialize the values.
xk = x0;
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
        lambdacur = lambda;
        xk = xk + lambdacur*(Afun(rxk,'transp'));
    elseif casel == 2
        % Landweber using line search.
        ATrk = Afun(rxk,'transp');
        lambda = norm(rxk)^2/norm(ATrk)^2;
        lambdacur = lambda;
        xk = xk + lambda*ATrk;
    elseif casel == 3
        % Landweber using psi1 or psi2.
        lambdacur = lambda(k);
        xk = xk + lambdacur*(Afun(rxk,'transp'));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rxk = b - Afun(xk,'notransp');
    
    % Check stopping rules:
    [stop,info,rk,dk] = check_stoprules(...
        stoprule, rxk, lambdacur, taudelta, k, kmax, rk, dk, dims);
        
    % If the current iteration is requested saved.
    if k == Knew(l+1) || stop
        l = l + 1;
        % Save the current iteration.
        % PCH rewrote the lines below.
        if strncmpi(stoprule,'ME',2) && stop && info(1)==3
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