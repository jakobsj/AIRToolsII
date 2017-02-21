function [A,b,x,theta,p,R,dw,sd] = ...
    fanbeamtomolinear(N,theta,p,R,dw,sd,isDisp)
%FANBEAMTOMOLINEAR Creates 2D fan-beam linear-detector tomography test problem.
%
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta,p)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta,p,R)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta,p,R,dw)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta,p,R,dw,sd)
%   [A,b,x,theta,p,R,d] = fanbeamtomolinear(N,theta,p,R,dw,sd,isDisp)
%
% This function creates a 2D tomography test problem with an N-times-N
% domain, using p rays in fan-formation for each angle in the vector theta.
% Unlike fanbeamtomo, in which the detector is curved, in fanbeamtomolinear
% the detector is linear, corresponding to the central slice of a
% flat-panel detector.
%
% Input:
%   N           Scalar denoting the number of discretization intervals in 
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the angles in degrees. Default: theta =
%               0:1:359.
%   p           Number of rays for each angle. Default: p =
%               round(sqrt(2)*N).
%   R           The distance from the source to the center of the domain
%               is R*N. Default: R = 2.
%   dw          Total width of linear detector. Same units as R.
%   sd          Source to detector distance. Same units as R.
%   isDisp      If isDisp is non-zero it specifies the time in seconds 
%               to pause in the display of the rays. If zero (the default), 
%               no display is shown.
%
% Output:
%   A           Coefficient matrix with N^2 columns and nA*p rows, 
%               where nA is the number of angles, i.e., length(theta).
%   b           Vector containing the rhs of the test problem.
%   x           Vector containing the exact solution, with elements
%               between 0 and 1.
%   theta       Vector containing the used angles in degrees.
%   p           The number of used rays for each angle.
%   R           The radius in side lengths. 
%   d           The span of the rays.
%
% See also: paralleltomo, fanbeamtomo, seismictomo.

% Modified to do linear detector array, from fanbeamtomo in AIRtools 1.0.
% Jakob Sauer Joergensen, 2012-04-04, DTU Compute, jakj@dtu.dk.

% Jakob Heide Joergensen, Maria Saxild-Hansen and Per Christian Hansen,
% June 21, 2011, DTU Informatics.

% Reference: A. C. Kak and M. Slaney, Principles of Computerized
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Default illustration:
if nargin < 7 || isempty(isDisp)
    isDisp = 0;
end

% Default value of sd.
if nargin < 6 || isempty(sd)
    sd = 3;
end
sd = sd*N;

% Default value of dw.
if nargin < 5 || isempty(dw)
    dw = 2.5;
end
dw = dw*N;

% Default value of R.
if nargin < 4 || isempty(R)
    R = 2;
end
R = R*N;

% Default value of the number of rays p.
if nargin < 3 || isempty(p)
    p = round(sqrt(2)*N);
end

% Default value of the angles theta.
if nargin < 2 || isempty(theta)
    theta = 0:359;
end

% Input check. The source must lie outside the domain.
if R < sqrt(2)/2*N
    error('R must be greater than half squareroot 2')
end

% Make sure theta is double precison to prevent round-off issues caused by
% single input.
theta = double(theta);

% Anonymous function rotation matrix
Omega_x = @(omega_par) [cosd(omega_par) -sind(omega_par)];
Omega_y = @(omega_par) [sind(omega_par)  cosd(omega_par)];

% nA denotes the number of angles.
nA = length(theta);

% The starting values both the x and the y coordinates.
x0 = 0;
y0 = R;
xy0 = [x0;y0];

% Width of detector element
dew = dw/p;

% Set angles to match linear detector
if mod(p,2)
    depos = (1:(p-1)/2)*dew;
    deposall = [-depos(end:-1:1),0,depos];
    omega = 90-atand(sd./depos);
    omega = [-omega(end:-1:1),0,omega];
else
    depos = ((1:p/2)-0.5)*dew;
    deposall = [-depos(end:-1:1),depos];
    omega = 90-atand(sd./depos);
    omega = [-omega(end:-1:1),omega];
end

% The intersection lines.
x = (-N/2:N/2)';
y = x;

% Initialize vectors that contains the row numbers, the column numbers and
% the values for creating the matrix A effiecently.
rows = zeros(2*N*nA*p,1);
cols = rows;
vals = rows;
idxend = 0;

% Prepare for illustration
if isDisp
    AA = rand(N);
    figure
end

% Loop over the chosen angles of the source.
for i = 1:nA
    
    % The starting points for the current angle theta.
    x0theta = Omega_x(theta(i))*xy0;
    y0theta = Omega_y(theta(i))*xy0;
    
    % The starting (center) direction vector (opposite xytheta) and
    % normalized.
    xytheta = [x0theta; y0theta];
    abtheta = -xytheta/R;
    
    % Illustration of the domain
    if isDisp % illustration of source
        clf
        pause(isDisp)
        imagesc((-N/2+.5):(N/2-0.5),(-N/2+.5):(N/2-0.5),AA), colormap gray,
        axis xy
        hold on
        axis(1.1*R*[-1,1,-1,1])
        axis equal
        plot(x0theta,y0theta,'o','color',[60 179 113]/255,...
            'linewidth',1.5,'markersize',10)
    end
    
    % Loop over the rays.
    for j = 1:p
        
        % The direction vector for the current ray.
        a = Omega_x(omega(j))*abtheta;
        b = Omega_y(omega(j))*abtheta;
        
        xdtheta = Omega_x(theta(i))*[deposall(j); R-sd];
        ydtheta = Omega_y(theta(i))*[deposall(j); R-sd];
        
        % Illustration of rays
        if isDisp
            plot([x0theta,xdtheta],[y0theta,ydtheta],'-',...
                'color',[220 0 0]/255,'linewidth',1.5)
            axis(R*[-1,1,-1,1])
        end
        
        % Use the parametrisation of line to get the y-coordinates of
        % intersections with x = k, i.e. x constant.
        tx = (x - x0theta)/a;
        yx = b*tx + y0theta;
        
        % Use the parametrisation of line to get the x-coordinates of
        % intersections with y = k, i.e. y constant.
        ty = (y - y0theta)/b;
        xy = a*ty + x0theta;
        
        % Collect the intersection times and coordinates.
        t = [tx; ty];
        xxy = [x; xy];
        yxy = [yx; y];
        
        % Sort the coordinates according to intersection time.
        [t, I] = sort(t);
        xxy = xxy(I);
        yxy = yxy(I);
        
        
        % Skip the points outside the box.
        I = (xxy >= -N/2 & xxy <= N/2 & yxy >= -N/2 & yxy <= N/2);
        xxy = xxy(I);
        yxy = yxy(I);
        
        % Skip double points.
        I = (abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
        xxy(I) = [];
        yxy(I) = [];
        
        % Illustration of the rays
        if isDisp
            set(gca,'Xticklabel',{})
            set(gca,'Yticklabel',{})
            pause(isDisp)
        end
        
        % Calculate the length within cell and determines the number of
        % cells which is hit.
        d = sqrt(diff(xxy).^2 + diff(yxy).^2);
        numvals = numel(d);
        
        % Store the values inside the box.
        if numvals > 0
            
            % Calculates the midpoints of the line within the cells.
            xm = 0.5*(xxy(1:end-1)+xxy(2:end));
            ym = 0.5*(yxy(1:end-1)+yxy(2:end));
            
            % Translate the midpoint coordinates to index.
            col = (floor(xm)+N/2)*N + (N - (floor(ym)+N/2));
            
            % Create the indices to store the values to vector for
            % later creation of A matrix.
            idxstart = idxend + 1;
            idxend = idxstart + numvals - 1;
            idx = idxstart:idxend;
            
            % Store row numbers, column numbers and values.
            rows(idx) = (i-1)*p + j;
            cols(idx) = col;
            vals(idx) = d;

        end
        
    end
    
end

% Truncate excess zeros.
rows = rows(1:idxend);
cols = cols(1:idxend);
vals = vals(1:idxend);

% Create sparse matrix A from the stored values.
A = sparse(rows,cols,vals,p*nA,N^2);

if nargout > 1
    
    % Create phantom head as a reshaped vector.
    x = myphantom(N);
    x = x(:);
    
    % Create rhs.
    b = A*x;
end

if nargout > 5
    R = R/N;
end

function X = myphantom(N)
%MYPHANTOM creates the modified Shepp-Logan phantom
%   X = myphantom(N)
% 
% This function create the modifed Shepp-Logan phantom with the
% discretization N x N, and returns it as a vector.
%
% Input:
%   N    Scalar denoting the nubmer of discretization intervals in each
%        dimesion, such that the phantom head consists of N^2 cells.
% 
% Output:
%   X    The modified phantom head reshaped as a vector

% This head phantom is the same as the Shepp-Logan except the intensities
% are changed to yield higher contrast in the image.
%
% Peter Toft, "The Radon Transform - Theory and Implementation", PhD
% thesis, DTU Informatics, Technical University of Denmark, June 1996.

%         A    a     b    x0    y0    phi
%        ---------------------------------
e =    [  1   .69   .92    0     0     0   
        -.8  .6624 .8740   0  -.0184   0
        -.2  .1100 .3100  .22    0    -18
        -.2  .1600 .4100 -.22    0     18
         .1  .2100 .2500   0    .35    0
         .1  .0460 .0460   0    .1     0
         .1  .0460 .0460   0   -.1     0
         .1  .0460 .0230 -.08  -.605   0 
         .1  .0230 .0230   0   -.606   0
         .1  .0230 .0460  .06  -.605   0   ];

xn = ((0:N-1)-(N-1)/2)/((N-1)/2);
Xn = repmat(xn,N,1);
Yn = rot90(Xn);
X = zeros(N);
     
% For each ellipse to be added     
for i = 1:size(e,1)
    a2 = e(i,2)^2;
    b2 = e(i,3)^2;
    x0 = e(i,4);
    y0 = e(i,5);
    phi = e(i,6)*pi/180;
    A = e(i,1);
    
    x = Xn-x0;
    y = Yn-y0;
    
    index = find(((x.*cos(phi) + y.*sin(phi)).^2)./a2 + ...
        ((y.*cos(phi) - x.*sin(phi))).^2./b2 <= 1);

    % Add the amplitude of the ellipse
    X(index) = X(index) + A;

end

% Return as vector and ensure nonnegative elements.
X = X(:); X(X<0) = 0;