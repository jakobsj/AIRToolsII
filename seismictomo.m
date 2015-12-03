function [A,b,x,s,p] = seismictomo(N,s,p,isDisp)
%SEISMICTOMO Creates a 2D seismic travel-time tomography test problem
%
%   [A,b,x,s,p] = seismictomo(N)
%   [A,b,x,s,p] = seismictomo(N,s)
%   [A,b,x,s,p] = seismictomo(N,s,p)
%   [A,b,x,s,p] = seismictomo(N,s,p,isDisp)
%
% This function creates a 2D seismic tomography test problem with an
% N-times-N domain, using s sources located on the right boundary and p
% receivers (seismographs) scattered along the left and top boundary.
% The rays are transmitted from each source to each receiver. 
%
% Inupt: 
%   N        Scalar denoting the number of discretization intervals in 
%            each dimesion, such that the domain consists of N^2 cells.
%   s        The number of sources in the right side of the domain
%            (default s = N).
%   p        The number of receivers (seismographs) equally spaced on the
%            surface and on the left side of the domain (default p = 2*N).
%   isDisp   If isDisp is nonzero it specifies the time in seconds 
%            to pause in the display of the rays. If zero (the default), 
%            no display is shown.
%
% Output:
%   A        Coefficient matrix with N^2 columns and s*p rows.
%   b        Vector containing the rhs of the test problem.
%   x        Vector containing the exact solution, with elements
%            between 0 and 1.
%   s        The number of sources.
%   p        The number of receivers (seismographs).
%
% See also: fanbeamtomo, paralleltomo.

% Jakob Sauer Jørgensen, Maria Saxild-Hansen and Per Chr. Hansen,
% October 1, 2014, DTU Compute.

if nargin < 4 || isempty(isDisp)
    isDisp = 0;
end

% Default number of sources.
if nargin < 2 || isempty(s)
    s = N;
end
if length(s)>1, error('s must be a scaler'), end

% Default number of receivers (seismographs).
if nargin < 3 || isempty(p)
    p = 2*N;
end

N2 = N/2;

% Determine the positions of the sources.
Ns = (N/s)/2;
x0 = N2*ones(s,1);                 
y0 = linspace(-N2+Ns,N2-Ns,s)';

% The intersection lines.
x = (-N2:N2)';
y = (-N2:N2)';

% The positions of the seismographs.
p1 = ceil(p/2); 	
p2 = floor(p/2); 
Np1 = (N/p1)/2;                 
Np2 = (N/p2)/2;

xp = [-N2*ones(p2,1); linspace(-N2+Np1,N2-Np1,p1)'];
yp = [linspace(-N2+Np2,N2-Np2,p2)'; N2*ones(p1,1)];

% Initialize the index for storing the results and the number of angles.
idxend = 0;
rows = zeros(2*N*s*p,1);
cols = rows;
vals = rows;

% Prepare for illustration.
if isDisp
    AA = rand(N);
    figure
end

% Loop over all the sources.
for i = 1:s
        
    % Illustation of the domain
    if isDisp

        clf
        pause(isDisp)
        imagesc((-N2+.5):(N2-0.5),(-N2+.5):(N2-0.5),AA), colormap gray,
        hold on
        axis xy
        axis equal
        axis([-N2 N2 -N2 N2])
        plot(x0,y0,'o','color',[60 179 113]/255,'linewidth',2,...
            'markersize',10)
        plot(xp,yp,'sb','MarkerSize',10,'linewidth',2)
        
    end
    
    % Determine the direction vector (a b).
    a = xp - x0(i);
    b = yp - y0(i);    
    
    % Loop over the seismographs.
    for j = 1:p
        
        % Use the parametrization of a line to get the y-coordinates of
        % intersections with x = constant.
        tx = (x - x0(i))/a(j);
        yx = b(j)*tx + y0(i);
        
        % Use the parametrisation of a line to get the y-coordinates of
        % intersections with y = constant.
        ty = (y - y0(i))/b(j);
        xy = a(j)*ty + x0(i);
        
        % Illustration of the rays.
        if isDisp
            plot(x,yx,'-','color',[220 0 0]/255,'linewidth',1.5)
            plot(xy,y,'-','color',[220 0 0]/255,'linewidth',1.5)
                                             
            set(gca,'Xticklabel',{})
            set(gca,'Yticklabel',{})
            pause(isDisp)
        end
        
        % Collect the intersection coordinates and times.
        t = [tx; ty];
        xxy = [x; xy];
        yxy = [yx; y];
        
        % Sort the times:
        [~,I] = sort(t);
        xxy = xxy(I);
        yxy = yxy(I);
        
        % Skip the points outside the box.
        I = find(xxy >= -N2 & xxy <= N2 & yxy >= -N2 & yxy <= N2);
        xxy = xxy(I);
        yxy = yxy(I);
        
        % Skip double points.
        I = find(abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
        xxy(I) = [];
        yxy(I) = [];
        
        % Calculate the lengths within the cells and the determine the
        % number of cells which are hit.
        d = sqrt(diff(xxy).^2 + diff(yxy).^2);
        numvals = numel(d);
        
        % Store the values inside the box.
        if numvals > 0
                           
            % Calculates the midpoints of the line within the cells.
            xm = 0.5*(xxy(1:end-1) + xxy(2:end)) + N2;
            ym = 0.5*(yxy(1:end-1) + yxy(2:end)) + N2;
            
            % Translate the midpoint coordinates to index.
            col = floor(xm)*N + (N - floor(ym));
            
            % Create the indices to store the values to a vector for
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
A = sparse(rows,cols,vals,s*p,N^2);

% Create the phantom.
if nargout > 1
    x = tectonic(N);
    b = A*x;
end

function x = tectonic(N)
% Creates a tectonic phantom of size N x N.

x = zeros(N);

N5 = round(N/5);
N13 = round(N/13);
N7 = round(N/7);
N20 = round(N/20);

% The right plate.
x(N5:N5+N7,5*N13:end) = 0.75;

% The angle of the right plate.
i = N5;
for j = 1:N20
    if rem(j,2) ~= 0
        i = i - 1;
        x(i,5*N13+j:end) = 0.75;
    end
end

% The left plate before the break.
x(N5:N5+N5,1:5*N13) = 1;

% The break from the left plate.
vector = N5:N5+N5;
for j = 5*N13:min(12*N13,N)
    if rem(j,2) ~= 0
        vector = vector + 1;
    end
    x(vector,j) = 1;
end

% Reshape the matrix to a vector.
x = x(:);