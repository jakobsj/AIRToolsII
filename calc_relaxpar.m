function [relaxpar, casel, sigma1tilde] = calc_relaxpar(relaxparinput, s1, kmax, atma, n)
%CALC_RELAXPAR Aux. function to compute relaxation parameter from input.
%
%   relaxpar = calc_relaxpar(relaxparinput)
%   [relaxpar, casel, sigma1tilde] = calc_relaxpar(relaxparinput, s1, kmax, atma, n)
%
% Short form is used by art and cart and sets the relaxpar to be used
% either to default values of 1 or 0.25, respectively, or assigns a value
% given by the user. If the value given by the user is outside the
% required interval a warning is given.
% 
% Long form is used by sirt and inputs and outputs are explained below.
% 
% Input:
%    relaxparinput     Any relaxation parameter or flag determining method 
%                      to use for determinin relaxpar as specified by user.
%    s1                The largest singular value if given by user.
%    kmax              The maximum number of SIRT iterations to run.
%    atma              A function handle to the iteration matrix
%                      characterizing the SIRT method, for which to compute
%                      the largest singular value.
%    n                 The number of columns in A.
%    
% Output:
%    relaxpar          The computed relaxation parameter or vector of
%                      iteration-dependent relaxation parameters.
%    casel             A flag indicating whether a constant lambda is
%                      returned (casel=1), line search is to be used
%                      (casel=2) or the psi1/psi2 strategies to be used
%                      (casel=3).
%    sigma1tilde       The computed largest singular value.
%
% See also: art.m, cart.m, sirt.m

% Jakob Sauer Jorgensen, Per Christian Hansen, Maria Saxild-Hansen
% 2017-03-03 DTU Compute

% First determine whether called from art, cart or sirt. Just testing for
% the specific method name such as cimmino would not cover the case of
% custom methods.
stack = dbstack;

switch stack(2).name
    
    case 'art'
        % Default choice 1. If user gave as input, use that value and 
        % throw warning if outside [0,2], but proceed.
        if isnan(relaxparinput)
            relaxpar = 1;
        else
            if relaxparinput <= 0 || relaxparinput >= 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The relaxpar value is outside the interval (0,2)');
            end
            relaxpar = relaxparinput;
        end
        
    case 'cart'
        % Default choice 1. If user gave as input, use that value and throw
        % warning if outside [0,2], but proceed.
        if isnan(relaxparinput)
            relaxpar = 0.25;
        else
            if relaxparinput <= 0 || relaxparinput >= 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The relaxpar value is outside the interval (0,2)');
            end
            relaxpar = relaxparinput;
        end
        
    case 'sirt'
        % Check if the largest singular value is given.
        if isnan(s1)
            % If not, calculate the largest singular value.
            optionsEIGS.disp = 0;
            sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
        else
            % Otherwise, use the given value.
            sigma1tilde = s1;
        end
        
        % Determine the relaxation parameter relaxpar.
        % If relaxparinput is nan, set default.
        if isnan(relaxparinput)
            
            % Define a default constant relaxpar value.
            relaxpar = 1.9/sigma1tilde^2;
            casel = 1;
            
            % If relaxpar is a scalar.
        elseif ~ischar(relaxparinput)
            
            % Checks if the given constant lambde value is unstable.
            if relaxparinput <= 0 || relaxparinput >= 2/sigma1tilde^2
                warning('MATLAB:UnstableRelaxParam',['The relaxpar value '...
                    'is outside the interval (0,%f)'],2/sigma1tilde^2)
            end
            relaxpar = relaxparinput;
            casel = 1;
            
        else
            % Calculate the relaxpar value according to the chosen method.
            if strncmpi(relaxparinput,'line',4)
                % Method: Line search
                casel = 2;
                
            elseif strncmpi(relaxparinput,'psi1',4)
                % Method: ENH psi1.
                casel = 3;
                
                % Precalculate the roots.
                z = calczeta(2:kmax-1);
                
                % Define the values for relaxpar according to the psi1
                % strategy modified or not.
                if strncmpi(relaxparinput,'psi1mod',7)
                    nu = 2;
                    relaxpar = [sqrt(2); sqrt(2); nu*2*(1-z)]/sigma1tilde^2;
                else
                    relaxpar = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tilde^2;
                end
                
            elseif strncmpi(relaxparinput,'psi2',4)
                % Method: ENH psi2.
                casel = 3;
                
                % Precalculate the roots.
                kk = 2:kmax-1;
                z = calczeta(kk);
                
                % Define the values for relaxpar according to the psi2
                % strategy modified or not.
                if strncmpi(relaxparinput,'psi2Mod',7)
                    nu = 1.5;
                    relaxpar = [sqrt(2); sqrt(2);
                        nu*2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
                else
                    relaxpar = [sqrt(2); sqrt(2);
                        2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
                end
                
            else
                error(['The chosen relaxation strategy is not valid '...
                    'for this method.'])
            end % end check of the class of relaxpar.
            
        end % end check of relaxpar strategies.
        
    otherwise
        error('This aux. function is only to be called from within art, cart or sirt.')
end

function z = calczeta(k)
%CALCZETA Calculates a specific root of a certain polynomial
%
%   z = calczeta(k)
%
% Calculates the unique root in the interval [0 1] of the polynomial equation
%       (2k -1)z^{k-1} - ( z^{k-2} + ... + z + 1 ) = 0
% by means of Newton's method.
%
% Input:
%   k  either a sorted vector or a scalar that determines the degree 
%      of the polynomial g(z).
% Output:
%   z  the found root(s).
%  
% See also: cav, cimmino, drop, landweber, sart, symkaczmarz.

% Maria Saxild-Hansen and Per Chr. Hansen, May 4, 2010, DTU Compute.

% Reference: T. Elfving, T. Nikazad, and P. C. Hansen, Semi-convergence and
% relaxation parameters for a class of SIRT algorithms, Elec. Trans. on
% Numer. Anal., 37 (201), pp. 321-336.

if k(1) < 2
    error('The polynomial is not defined for k < 2')
end

% The starting guess for Newton's method
z0 = 1;
% The number of Newton iterations
kmax = 6;
z = zeros(length(k),1);

% Finds the root with Newton's method for all the values in k
for i = 1:length(k)
    z(i) = myNewton(z0,kmax,k(i));
end

function [x,k] = myNewton(x0,kmax,grad)
% The Newton method specially defined for this problem

% Initialize the parameters
k = 0;
x = x0;
[f, df] = fung(x,grad);
h = f/df;

% Iterate using Newtoms method
while  k < kmax
    k = k+1;
    x = x - h;
    [f, df] = fung(x,grad);
    h = f/df;
end

function [f df] = fung(x,k)
% Uses Horner's algorithm to evaluate the ploynomial and its derivative

C = [0 2*k-1 -ones(1,k-1)];
f = C(1);
df = 0;

for i = 2:length(C)
    df = df*x + f;
    f = f*x + C(i);
end
