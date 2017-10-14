function tau = train_dpme(A,b_exact,x_exact,method,type,delta,s,kmax,options)
%TRAIN_DPME  Training method for the stopping rules DP and ME
%
%   tau = train_dpme(A,b_exact,x_exact,method,type,delta,s,kmax)
%   tau = train_dpme(A,b_exact,x_exact,method,type,delta,s,kmax,options)
%
% This function determines the parameter tau for a given method, for the
% stopping rule DP or ME, by training on a problem given by A, b_exact, and
% x_exact. The noise level is delta, s samples of the noisy right-hand side
% are created, and from these samples the value of tau is determined.
%
% Input:
%   A           m times n matrix or function handle to matrix-free version.
%   b_exact     m times 1 vector containing the exact rhs.
%   x_exact     n times 1 vectir containing the exact solution.
%   method      Function handle to a SIRT, ART or CART method.
%   type        String that should be either 'DP' for Discrepancy Principle
%               or 'ME' for the Monotone Error rule.  ME can only be chosen
%               if method is a SIRT-method.
%   delta       Scalar denoting the noise level || e ||_2.
%   s           The number of used samples.
%   kmax        Maximum number of iterations.
%   options     Struct used in the call of the method.  For this strategy
%               the field stoprule cannot be used.
%
% Output:
%   tau         Scalar containing the trained parameter.
%
% See also: demo_training, train_relaxpar.

% Reference: T. Elfving and T. Nikazad, Stopping rules for Landweber-type
% iteration, Inverse Problems, 23 (2007), pp. 1417-1432.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Remove any stopping rule given and make sure is set to none.
if nargin == 9
    if isfield(options,'stoprule')
        options = rmfield(options,'stoprule');
    end
end
options.stoprule.type = 'none';

% If A is not a function (i.e., A is a matrix or an object), convert A
% to a function.
if isa(A,'function_handle')
    Afun = A;
else
    Afun = @(xx,transp_flag) afun_matrix(xx,transp_flag,A);
end
m = Afun([],'size');
m = m(1);

% Generate s noisy samples of the rhs.
e = randn(m,s);
e = delta*e./(ones(m,1)*sqrt(sum(e.*e)));
b = b_exact*ones(1,s) + e;

% Initialize the tauh vector.
tauh = zeros(s,1);

% The maximum number of iterations.
if kmax < 10
    kstep = kmax;
else
    kstep = 10;
end

% Loop over all generated samples of the rhs b.
for i = 1:s
    
    % Prepare for solving with the i'th noisy rhs.
    bi = b(:,i);
    X = zeros(size(x_exact,1),1);
    x0 = zeros(size(x_exact));
    k = kstep;
    K = 1:k;
    Ek = [];
    
    stop = 0;
    
    % Iterate as long as the minimum residual is not found.
    while ~stop
        
        % Perform K iterations.
        Xny = method(Afun,bi,K,x0,options);
        
        % Compute the error and find the minimum.
        xd = Xny - repmat(x_exact,1,size(Xny,2));
        Ek = [Ek sqrt(sum(xd.*xd,1))];
        X = [X Xny];
        [~, I] = min(Ek);
        
        if I ~= k || k >= kmax
            % Stop the loop and define the optimal number of iterations.
            stop = 1;
            kopt = I;
        else
            % Continue iterations, update iterations and starting value.
            if k+kstep > kmax
                K = 1:(kmax-k);
                k = kmax;
            else
                k = k + kstep;
            end
            x0 = Xny(:,end);
        end
        
    end
    
    % Compute the residual for the optimal iteration and the previous one.
    rk = bi-Afun(X(:,kopt),'notransp');
    
    % Split in the types DP and ME.
    if strncmpi(type,'DP',2)
        tauk = norm(rk)/delta;
    elseif strncmpi(type,'ME',2)
        rkm1 = bi-Afun(X(:,kopt-1),'notransp');
        tauk = 0.5*(rkm1'*(rkm1+rk))/(delta*norm(rkm1));
    else
        error('The chosen type is invalid')
    end
    
    % Define the tau value for the current right-hand side.
    tauh(i) = tauk;
end

% Return the optimal tau value.
tau = mean(tauh);