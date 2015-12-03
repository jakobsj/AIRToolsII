function lambda = trainLambdaSIRT(A,b,x_ex,method,kmax,options)
%TRAINLAMBDASIRT Training to determine optimal lambda for SIRT method
%
%   lambda = trainLambdaSIRT(A,b,x_ex,method)
%   lambda = trainLambdaSIRT(A,b,x_ex,method,kmax)
%   lambda = trainLambdaSIRT(A,b,x_ex,method,kmax,options)
%
% This function determines the optimal value of the relaxation parameter
% lambda for one of the SIRT methods cav, cimmino, drop, landweber, and sart.
% The optimal value of lambda is defined as the value that gives rise to
% the fastest convergence to the smallest relative error in the solution. 
%
% Input:
%   A           m times n matrix.
%   b           m times 1 vector containing the right-hand side.
%   x_ex        n times 1 vector containing the exact solution.
%   method      Function handle to one of the SIRT methods.
%   kmax        Scalar that determines the maximum number of iterations
%               of the used SIRT-method. Default value is 1000.
%   options     Struct used in the call of the method. For this strategy
%               the fields stoprule and lambda cannot be used.
%
% Output:   
%   lambda      Scalar containing the found lambda value.
%
% See also: trainLambdaART

% Maria Saxild-Hansen and Per Chr. Hansen, June 10, 2010, DTU Compute.

% Input check for options. Ensure that no stoprule type is chosen.
if nargin == 6
    options.stoprule.type = 'none';
    if isfield(options,'lambda')
        options = rmfield(options,'lambda');
    end
else
    options = [];
end 

% Input checks and default values:
if nargin < 5 || isempty(kmax)
    % Default maximum number of iterations.
    kmax = 1000;
end

% Minimum error interval.
pct = 0.01;

% Initialize the number of iterations.
stepk = 10;
K = 1:stepk;
k = 0;
stop = 0;
rEr = [];
x0 = [];

% Define the norm of the exact solution.
normxex = norm(x_ex);

% Step 1: Determine the minimum error at the default value.
while ~stop
    % Ensure only the maximum number of iterations are performed.
    if k+stepk > kmax
        K = 1:(kmax-k);
        k = kmax;
    else
        % Update the number of iterations.
        k = k + stepk;
    end
        
    % If the first run in the loop.
    if k <= stepk
        
        switch func2str(method)
            case {'landweber','cimmino','cav','drop'}
                
                % Determine the solutions to the K values.
                [Xnew info restart] = method(A,b,K,x0,options);
                
                % Assign the restart struct to options, such that M, T and
                % s1 are not recalculated in the next call of the method.
                options.restart = restart;
                lambdamax = 2/restart.s1^2;
                
            case {'sart'}
                
                % Determine the solutions to the K values.
                [Xnew info restart] = method(A,b,K,x0,options);
                
                % Assign the restart struct to options, such that M, T and
                % s1 are not recalculated in the next call of the method.
                options.restart = restart;
                lambdamax = 2;
                
            otherwise
                error(['Unknown method ',func2str(method)])
        end

    else
        
        % Determine the solutions to the next K values.
        Xnew = method(A,b,K,x0,options);
    end
    
    % Calculate the relative error for the new solutions and collect all
    % the relative errors.
    deltax = Xnew - repmat(x_ex,1,size(Xnew,2));
    rEr = [rEr sqrt(sum(deltax.*deltax,1))/normxex];
    
    % Find the minimum relative error.
    [minE kopt] = min(rEr);
    
    % If the minimum relative error is in the last iteration, then run the
    % loop again with stepk new iterations.  Else the minimum error is
    % found and the loop is stopped.
    if kopt == k && k < kmax
        x0 = Xnew(:,end);
    else
        stop = 1;
        
    end
end

% Determine the minimum error interal.
pinter = [minE-pct*minE minE+pct*minE];

% Step 2: Find the lambda for which the minimum relative error is
% within pinter and where kopt is the smallest.

% Define a lambda value near the end of the allowed interval, calculate
% the relative errors for the solutions, and find the minimum.
options.lambda = lambdamax-0.01*lambdamax;
Xend = method(A,b,1:kopt,[],options);
deltaxend = Xend - repmat(x_ex,1,kopt);
rErend = sqrt(sum(deltaxend.*deltaxend,1))/normxex;

[minEend koptend] = min(rErend);

% If kopt is equal to koptend, then the global minimum error is not reached,
% hence all the k values are equal.  The best lambda value is determined by
% the minimum relative error.

% Define the search interval.
r = (3-sqrt(5))/2;

i1 = 0;
i2 = lambdamax;
i3 = i1 + r*(i2-i1);
i4 = i1 + (1-r)*(i2-i1);

% Determine the function values for i3 and i4.
options.lambda = i3;
xnew = method(A,b,1:kopt,[],options);
if kopt == koptend
    f3 = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
    
    % Define x3 always to be inside the interval.
    x3 = pinter(2);
else
    [x3 f3] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
    
end

options.lambda = i4;
xnew = method(A,b,1:kopt,[],options);
if kopt == koptend
    f4 = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
    
    % Define x4 always to be inside the interval.
    x4 = pinter(2);
else
    [x4 f4] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
end

while abs(i3-i4) > lambdamax*0.01
    
    % Ensure that the minimum error is reached in the left part of
    % the interval.
    if x3 > pinter(2)
        
        % i1 is removed from the search interval.
        z = i3 + (1-r)*(abs(i2-i3));
        
        % [i1 i3 i4 i2] <---- [i3 i4 z i2]
        i1 = i3;
        i3 = i4;
        f3 = f4;
        x3 = x4;
        i4 = z;
        
        % Find the new function value for f4.
        options.lambda = i4;
        xnew = method(A,b,1:kopt,[],options);
        [x4 f4] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        
    elseif x4 > pinter(2)
        % Ensure that the minimum error is reached in the left part of
        % the interval.
        
        % i2 is removed from the search interval.
        z = i1 + r*(abs(i4-i1));
        
        % [i1 i3 i4 i2] <---- [i1 z i3 i4]
        i2 = i4;
        i4 = i3;
        f4 = f3;
        
        x4 = x3;
        i3 = z;
        
        % Find the new value for f3.
        options.lambda = i3;
        xnew = method(A,b,1:kopt,[],options);
        [x3 f3] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        
    elseif f3 >= f4
        
        % If the relative error for f3 is greater than or equal to f4, then
        % i1 is removed from the interval.
        z = i3 + (1-r)*(abs(i2-i3));
        
        % [i1 i3 i4 i2] <---- [i3 i4 z i2]
        i1 = i3;
        i3 = i4;
        f3 = f4;
        x3 = x4;
        i4 = z;
        
        % Find the new value for f4
        options.lambda = i4;
        xnew = method(A,b,1:kopt,[],options);
        if kopt == koptend
            f4 = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        else
            [x4 f4] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        end
        
    else
        
        % i2 is removed from the search interval.
        z = i1 + r*(abs(i4-i1));
        
        % [i1 i3 i4 i2] <---- [i1 z i3 i4]
        i2 = i4;
        i4 = i3;
        f4 = f3;
        x4 = x3;
        i3 = z;
        
        % Find the new value for f3.
        options.lambda = i3;
        xnew = method(A,b,1:kopt,[],options);
        if kopt == koptend
            f3 = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        else
            [x3 f3] = min(sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1))/normxex);
        end
        
    end
end

% Define the optimal lambda value.
lambda = (i3+i4)/2;