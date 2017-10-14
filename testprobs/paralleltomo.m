function [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp,isMatrix)
%PARALLELTOMO  Creates a 2D parallel-beam tomography test problem
%
%   [A,b,x,theta,p,d] = paralleltomo(N)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp)
%   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp,isMatrix)
%
% This function uses the "line model" to create a 2D X-ray tomography test
% problem with an N-times-N pixel domain, using p parallel rays for each
% angle in the vector theta.
%
% Input:
%   N         Scalar denoting the number of discretization intervals in
%             each dimesion, such that the domain consists of N^2 cells.
%   theta     Vector containing the projetion angles in degrees.
%             Default: theta = 0:1:179.
%   p         Number of rays for each angle. Default: p = round(sqrt(2)*N).
%   d         Scalar denoting the distance from the first ray to the last.
%             Default: d = p-1.
%   isDisp    If isDisp is non-zero it specifies the time in seconds
%             to pause in the display of the rays. If zero (the default),
%             no display is shown.
%   isMatrix  If non-zero, a sparse matrix is returned in A (default).
%             If zero, instead a function handle is returned.
%
% Output:
%   A         If input isMatrix is 1 (default): coefficient matrix with
%             N^2 columns and length(theta)*p rows.
%             If input isMatrix is 0: A function handle representing a
%             matrix-free version of A in which the forward and backward
%             operations are computed as A(x,'notransp') and A(y,'transp'),
%             respectively, for column vectors x and y of appropriate size.
%             The size of A can be retrieved using A([],'size'). The matrix
%             is never formed explicitly, thus saving memory.
%             elements in A.
%   b         Vector containing the rhs of the test problem.
%   x         Vector containing the exact solution, with elements
%             between 0 and 1.
%   theta     Vector containing the used angles in degrees.
%   p         The number of used rays for each angle.
%   d         The distance between the first and the last ray.
%
% See also: fancurvedtomo, fanlineartomo, seismictomo, seismicwavetomo,
%           sphericaltomo.

% Reference: A. C. Kak and M. Slaney, Principles of Computerized
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Default value of the projection angles theta.
if nargin < 2 || isempty(theta)
    theta = 0:179;
end

% Make sure theta is double precison to prevent round-off issues caused by
% single input.
theta = double(theta);

% Default value of the number of rays.
if nargin < 3 || isempty(p)
    p = round(sqrt(2)*N);
end

% Default value of d.
if nargin < 4 || isempty(d)
    d = p-1;
end

% Default illustration: False
if nargin < 5 || isempty(isDisp)
    isDisp = 0;
end

% Default matrix or matrix-free function handle? Matrix
if nargin < 6 || isempty(isMatrix)
    isMatrix = 1;
end

% Construct either matrix or function handle
if isMatrix
    A = get_or_apply_system_matrix(N,theta,p,d,isDisp);
else
    A = @(U,TRANSP_FLAG) get_or_apply_system_matrix(N,theta,p,d,isDisp,...
        U,TRANSP_FLAG);
end

if nargout > 1
    % Create phantom head as a reshaped vector.
    x = phantomgallery('shepplogan',N);
    x = x(:);
    % Create rhs.
    if isMatrix
        b = A*x;
    else
        b = A(x,'notransp');
    end
end


function A = get_or_apply_system_matrix(N,theta,p,d,isDisp,u,transp_flag)

% Define the number of angles.
nA = length(theta);

% The starting values both the x and the y coordinates.
x0 = linspace(-d/2,d/2,p)';
y0 = zeros(p,1);

% The intersection lines.
x = (-N/2:N/2)';
y = x;

% Prepare for illustration
if isDisp
    AA = phantomgallery('smooth',N);
    figure
end

% Deduce whether to set up matrix or apply to input, depending on whether
% input u is given.
isMatrix = (nargin < 6);

if isMatrix
    
    % Initialize vectors that contains the row numbers, the column numbers
    % and the values for creating the matrix A effiecently.
    rows = zeros(2*N*nA*p,1);
    cols = rows;
    vals = rows;
    idxend = 0;
    
    II = 1:nA;
    JJ = 1:p;
else
    % If transp_flag is 'size', only return size of operator, otherwise set
    % proper size of output.
    switch transp_flag
        case 'size'
            A = [p*nA,N^2];
            return
        case 'notransp' % Forward project.
            if length(u) ~= N^2
                error('Incorrect length of input vector')
            end
            A = zeros(p*nA,1);
        case 'transp' % Backproject
            if length(u) ~= p*nA
                error('Incorrect length of input vector')
            end
            A = zeros(N^2,1);
    end
    
    % If u is a Cartesian unit vector then we only want to compute a single
    % row of A; otherwise we want to multiply with A or A'.
    if strcmpi(transp_flag,'transp') && nnz(u) == 1 && sum(u) == 1
        % Want to compute a single row of A, stored as a column vector.
        ell = find(u==1);
        JJ = mod(ell,p);  if JJ==0, JJ = p; end
        II = (ell-JJ)/p + 1;
    else
        % Want to multiply with A or A'.
        II = 1:nA;
        JJ = 1:p;
    end
end

% Loop over the chosen angles.
for i = II
    
    % Illustration of the domain.
    if isDisp
        clf
        pause(isDisp)
        imagesc((-N/2+.5):(N/2-0.5),(-N/2+.5):(N/2-0.5),flipud(AA))
        colormap gray
        axis xy
        hold on
        axis equal
        axis(0.7*[-N N -N N])
    end
    
    % All the starting points for the current angle.
    x0theta = cosd(theta(i))*x0-sind(theta(i))*y0;
    y0theta = sind(theta(i))*x0+cosd(theta(i))*y0;
    
    % The direction vector for all rays corresponding to the current angle.
    a = -sind(theta(i));
    b = cosd(theta(i));
    
    % Loop over the rays.
    for j = JJ
        
        % Use the parametrisation of line to get the y-coordinates of
        % intersections with x = constant.
        tx = (x - x0theta(j))/a;
        yx = b*tx + y0theta(j);
        
        % Use the parametrisation of line to get the x-coordinates of
        % intersections with y = constant.
        ty = (y - y0theta(j))/b;
        xy = a*ty + x0theta(j);
        
        % Illustration of the rays.
        if isDisp
            
            plot(x,yx,'-','color',[220 0 0]/255,'linewidth',1.5)
            plot(xy,y,'-','color',[220 0 0]/255,'linewidth',1.5)
            
            set(gca,'Xtick',[],'Ytick',[])
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
        col = [];
        
        % Store the values inside the box.
        if numel(aval) > 0
            
            % If the ray is on the boundary of the box in the top or to the
            % right the ray does not by definition lie with in a valid cell.
            if ~((b == 0 && abs(y0theta(j) - N/2) < 1e-15) || ...
                    (a == 0 && abs(x0theta(j) - N/2) < 1e-15)       )
                
                % Calculates the midpoints of the line within the cells.
                xm = 0.5*(xxy(1:end-1)+xxy(2:end)) + N/2;
                ym = 0.5*(yxy(1:end-1)+yxy(2:end)) + N/2;
                
                % Translate the midpoint coordinates to index.
                col = floor(xm)*N + (N - floor(ym));
                    
            end
        end
        if ~isempty(col)
            if isMatrix
                % Create the indices to store the values to vector for
                % later creation of A matrix.
                idxstart = idxend + 1;
                idxend = idxstart + numel(col) - 1;
                idx = idxstart:idxend;
                
                % Store row numbers, column numbers and values.
                rows(idx) = (i-1)*p + j;
                cols(idx) = col;
                vals(idx) = aval;
            else
                % If any nonzero elements, apply forward or back operator
                
                if strcmp(transp_flag,'notransp')
                    % Insert inner product with u into w.
                    A( (i-1)*p+j ) = aval'*u(col);
                else  % Adjoint operator.
                    A(col) = A(col) + u( (i-1)*p+j )*aval;
                end
            end
        end
    end % end j
end % end i

if isMatrix
    % Truncate excess zeros.
    rows = rows(1:idxend);
    cols = cols(1:idxend);
    vals = vals(1:idxend);
    
    % Create sparse matrix A from the stored values.
    A = sparse(rows,cols,vals,p*nA,N^2);
end
