function w = OPparalleltomo(u,transp_flag,N,theta,p,d)
%OPparalleltomo  Operator corresponding to the test problem paralleltomo
%
%   w = OPparalleltomo(u,transp_flag,N)
%   w = OPparalleltomo(u,transp_flag,N,theta)
%   w = OPparalleltomo(u,transp_flag,N,theta,p)
%   w = OPparalleltomo(u,transp_flag,N,theta,p,d)
%
% This function applies forward operator to a vector, the adjoint/
% transposed to a vector, or returns the size of the operator. The operator
% corresponds to the system matrix A as set up by the paralleltomo function.
%
% Input:
%   u           If input transp_flag is 'notransp': image vector of length N^2.
%               If input transp_flag is 'transp': sinogram vector of size nA*p.
%               If input transp_flag is 'size': Unused.
%   transp_flag String specifying which operator to apply. If 'size', simply
%               return the dimensions of the operator; if 'notransp', apply
%               the forward operator (projection) to the input u; if 'transp', 
%               apply the adjoint (backprojection) operator to input u.
%   N           Scalar denoting the number of discretization intervals in 
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the angles in degrees. Default: theta = 
%               0:1:179.
%   p           Number of parallel rays for each angle. Default: p =
%               round(sqrt(2)*N).
%   d           Scalar denoting the distance from the first ray to the last.
%               Default: d = sqrt(2)*N.
%
% Output:
%   w           Either A*u, A'*u, or size(A).
% 
% See also: OPfanbeamtomo

% Jakob Sauer Joergensen, Maria Saxild-Hansen and Per Christian Hansen,
% November 7, 2015, DTU Compute.

% Reference: A. C. Kak and M. Slaney, Principles of Computerized 
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Default value of the angles theta.
if nargin < 4 || isempty(theta)
    theta = 0:179;
end

% Default value of the number of rays.
if nargin < 5 || isempty(p)
    p = round(sqrt(2)*N);
end

% Default value of d.
if nargin < 6 || isempty(d)
    d = p-1;
end

% Define the number of angles.
nA = length(theta);

% If transp_flag is 'size', only return size of operator, otherwise set
% proper size of output.
switch transp_flag
    case 'size'
        w = [p*nA,N^2];
        return
    case 'notransp' % Forward project.
        if length(u) ~= N^2, error('Incorrect length of input vector'), end
        w = zeros(p*nA,1);
    case 'transp' % Backproject
        if length(u) ~= p*nA, error('Incorrect length of input vector'), end
        w = zeros(N^2,1);
end

% The starting values both the x and the y coordinates. 
x0 = linspace(-d/2,d/2,p)';
y0 = zeros(p,1);

% The intersection lines.
x = (-N/2:N/2)';
y = x;

if strcmpi(transp_flag,'transp') && nnz(u) == 1 && sum(u) == 1
    % Compute a single row of A, stored as a column vector.
    ell = find(u==1);
    J = mod(ell,p);  if J==0, J = p; end
    I = (ell-J)/p + 1;
else
    I = 1:nA;
    J = 1:p;
end

% Loop over the chosen angles.
for i = I    
        
    % All the starting points for the current angle.
    x0theta = cosd(theta(i))*x0-sind(theta(i))*y0;
    y0theta = sind(theta(i))*x0+cosd(theta(i))*y0;
    
    % The direction vector for all the rays corresponding to the current 
    % angle.
    a = -sind(theta(i));
    b = cosd(theta(i));
    
    % Loop over the rays.
    for j = J
        
        % Use the parametrisation of line to get the y-coordinates of
        % intersections with x = k, i.e. x constant.
        tx = (x - x0theta(j))/a;
        yx = b*tx + y0theta(j);
        
        % Use the parametrization of line to get the x-coordinates of
        % intersections with y = k, i.e. y constant.
        ty = (y - y0theta(j))/b;
        xy = a*ty + x0theta(j);
                
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
        v = sqrt(diff(xxy).^2 + diff(yxy).^2);
        numvals = numel(v);
        col = [];
        
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
                
            end
        end
        
        % If any nonzero elements, apply forward or back operator
        if ~isempty(col)
            if strcmp(transp_flag,'notransp')
                % Insert inner product with u into w.
                w( (i-1)*p+j ) = v'*u(col);
            else  % Adjoint operator.
                w(col) = w(col) + u( (i-1)*p+j )*v;
            end
        end
   end
              
end