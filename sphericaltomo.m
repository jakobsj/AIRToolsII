function [A,b,x,angles,numCircles] = sphericaltomo(N,angles,numCircles)
%STtomo Creates a 2D spherical Radon tomography test problem
%
%   [A,b,x,angles,numCircles] = sphericaltomo(N)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles,numCircles)
%
% This function genetates a tomography test problem based on the spherical
% Radon tranform where data consists of integrals along circles.  This type
% of problem arises, e.g., in photoacoustic imaging.
%
% The image is square whose center is at the origo.  The centers for the
% integration circles are placed on a circle just outside the square.
% For each circle center we integrate along a number of circles with
% different radii, using the periodic trapezoidal rule.
%
% Input:
%   N           Scalar denoting the number of pixels in each dimesion, such
%               that the image domain consists of N^2 cells.
%   angles      Vector containing the angles to the circle centers in degrees.
%               Default: angles = 0:2:358.
%   numCircles  Number of concentric integration circles for each center.
%               Default: numCircles = round(sqrt(2)*N).
%
% Output:
%   A           Sparse coefficient matrix with N^2 columns and
%               length(angles)*numColumns rows.
%   b           Vector containing the rhs of the test problem.
%   x           Vector containing the exact solution, with elements
%               between 0 and 1.
%   angles      Similar to input.
%   numCircles  Similar to input.
%
% See also: xxx

% Per Christian Hansen, DTU Compute, May 15, 2017.
% Inspired by Matlab code written by Jürgen Frikel.

% Default value of the angles to the circle centers.
if nargin < 2 || isempty(angles)
    angles = 0:2:358;
end

% Default value of the number of circles.
if nargin < 3 || isempty(numCircles)
    numCircles = round(sqrt(2)*N);
end

% Radii for the circles.
radii  = linspace(0,2,numCircles+1);
radii  = radii(2:end);

% Image coordinates.
centerImg = ceil(N/2);

% Determine the quarature parameters.
dx = sqrt(2)/N;
nPhi = ceil((4*pi/dx)*radii);
dPhi = 2*pi./nPhi;

ii = zeros(0,1);
jj = zeros(0,1);
aa = zeros(0,1);

% Loop over angles and circles.
for m=1:length(angles)
    xix = cosd(angles(m));
    xiy = sind(angles(m));
    for n=1:numCircles
            k = (0:nPhi(n))*dPhi(n);
            col = round( (xix + radii(n)*cos(k))/dx + centerImg );
            row = round( (xiy + radii(n)*sin(k))/dx + centerImg );
            II = find(col>0 & col<=N & row>0 & row<=N);
            row = row(II);
            col = col(II);
            J = row + (col-1)*N;
            Ju = unique(J);
            w = histc(J,Ju);
            i = m + (n-1)*length(angles);
            ii = [ii;repmat(i,length(Ju),1)];
            jj = [jj;Ju(:)];
            aa = [aa;(2*pi*radii(n)/nPhi(n))*w(:)];
    end
end
A = sparse(ii,jj,aa,length(angles)*numCircles,N^2);

if nargout > 1
    x = phantomgallery('shepplogan',N);
    x = x(:);
    b = A*x;
end