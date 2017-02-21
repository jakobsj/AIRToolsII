function [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp)
%PARALLELTOMO Creates a 2D tomography test problem using parallel beams
%
%   [A,b,x,theta,p,d] = paralleltomo(N)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp)
%
% This function creates a 2D tomography test problem with an N-times-N
% domain, using p parallel rays for each angle in the vector theta.
%
% Input: 
%   N           Scalar denoting the number of discretization intervals in 
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the projetion angles in degrees.
%               Default: theta = 0:1:179.
%   p           Number of parallel rays for each angle.
%               Default: p = round(sqrt(2)*N).
%   d           Scalar denoting the distance from the first ray to the last.
%               Default: d = p-1.
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
%   d           The distance between the first and the last ray.
% 
% See also: fanbeamtomo, seismictomo.

% Jakob Sauer J�rgensen, Maria Saxild-Hansen and Per Christian Hansen,
% Nov. 5, 2015, DTU Compute.

% Reference: A. C. Kak and M. Slaney, Principles of Computerized 
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Default value of the projection angles theta.
if nargin < 2 || isempty(theta)
    theta = 0:179;
end

% Default value of the number of rays.
if nargin < 3 || isempty(p)
    p = round(sqrt(2)*N);
end

% Default value of d.
if nargin < 4 || isempty(d)
    d = p-1;
end

% Default illustration:
if nargin < 5 || isempty(isDisp)
    isDisp = 0;
end

% Make sure theta is double precison to prevent round-off issues caused by
% single input.
theta = double(theta);

% Define the number of angles.
nA = length(theta);

% The starting values both the x and the y coordinates. 
x0 = linspace(-d/2,d/2,p)';
y0 = zeros(p,1);

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

% Loop over the chosen angles.
for i = 1:nA    
    
    % Illustration of the domain
    if isDisp
        clf
        pause(isDisp)
        imagesc((-N/2+.5):(N/2-0.5),(-N/2+.5):(N/2-0.5),AA), colormap gray,
        hold on
        axis xy
        axis equal
        axis([-N/2-1 N/2+1 -N/2-1 N/2+1])        
    end
    
    % All the starting points for the current angle.
    x0theta = cosd(theta(i))*x0-sind(theta(i))*y0;
    y0theta = sind(theta(i))*x0+cosd(theta(i))*y0;
    
    % The direction vector for all the rays corresponding to the current 
    % angle.
    a = -sind(theta(i));
    b = cosd(theta(i));
    
    % Loop over the rays.
    for j = 1:p
        
        % Use the parametrisation of line to get the y-coordinates of
        % intersections with x = constant.
        tx = (x - x0theta(j))/a;
        yx = b*tx + y0theta(j);
        
        % Use the parametrisation of line to get the x-coordinates of
        % intersections with y = constant.
        ty = (y - y0theta(j))/b;
        xy = a*ty + x0theta(j);

        % Illustration of the rays
        if isDisp           
            
            plot(x,yx,'-','color',[220 0 0]/255,'linewidth',1.5)
            plot(xy,y,'-','color',[220 0 0]/255,'linewidth',1.5)
                     
            set(gca,'Xticklabel',{})
            set(gca,'Yticklabel',{})
            pause(isDisp)
        end
        
        % Collect the intersection times and coordinates. 
        t = [tx; ty];
        xxy = [x; xy];
        yxy = [yx; y];
        
        % Sort the coordinates according to intersection time.
        [~,I] = sort(t);
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
        
        % Calculate the length within cell and determines the number of
        % cells which is hit.
        aval = sqrt(diff(xxy).^2 + diff(yxy).^2);
        numvals = numel(aval);
        
        % Store the values inside the box.
        if numvals > 0
            
            % If the ray is on the boundary of the box in the top or to the
            % right the ray does not by definition lie with in a valid cell.
            if ~((b == 0 && abs(y0theta(j) - N/2) < 1e-15) || ...
                 (a == 0 && abs(x0theta(j) - N/2) < 1e-15)       )
                
                % Calculates the midpoints of the line within the cells.
                xm = 0.5*(xxy(1:end-1)+xxy(2:end)) + N/2;
                ym = 0.5*(yxy(1:end-1)+yxy(2:end)) + N/2;
                
                % Translate the midpoint coordinates to index.
                col = floor(xm)*N + (N - floor(ym));
                
                % Create the indices to store the values to vector for
                % later creation of A matrix.
                idxstart = idxend + 1;
                idxend = idxstart + numvals - 1;
                idx = idxstart:idxend;
                
                % Store row numbers, column numbers and values. 
                rows(idx) = (i-1)*p + j;
                cols(idx) = col;
                vals(idx) = aval;   
                
            end
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
    x = phantomgallery('shepplogan',N);
    % Create rhs.
    b = A*x;
end