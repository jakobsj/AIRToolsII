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