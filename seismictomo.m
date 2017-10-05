function [A,b,x,s,p] = seismictomo(N,s,p,isDisp,isMatrix)
%SEISMICTOMO  Creates a 2D seismic travel-time tomography test problem
%
%   [A,b,x,s,p] = seismictomo(N)
%   [A,b,x,s,p] = seismictomo(N,s)
%   [A,b,x,s,p] = seismictomo(N,s,p)
%   [A,b,x,s,p] = seismictomo(N,s,p,isDisp)
%   [A,b,x,s,p] = seismictomo(N,s,p,isDisp,isMatrix)
%
% This function uses the "line model" to create a 2D seismic travel-time
% tomography test problem with an N-times-N pixel domain, using s sources
% located on the right boundary and p receivers (seismographs) scattered
% along the left and top boundary.  The rays are transmitted from each
% source to each receiver. 
%
% Inupt: 
%   N         Scalar denoting the number of discretization intervals in 
%             each dimesion, such that the domain consists of N^2 cells.
%   s         The number of sources in the right side of the domain.
%             Default s = N.
%   p         The number of receivers (seismographs) equally spaced on the
%             surface and on the left side of the domain. Default p = 2*N.
%   isDisp    If isDisp is nonzero it specifies the time in seconds 
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
%   b         Vector containing the rhs of the test problem.
%   x         Vector containing the exact solution, with elements
%             between 0 and 1.
%   s         The number of sources.
%   p         The number of receivers (seismographs).
%
% See also: paralleltomo, fancurvedtomo, fanlineartomo, seismicwavetomo,
%           sphericaltomo.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

if nargin < 5 || isempty(isMatrix)
    isMatrix = 1;
end

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

% Construct either matrix or function handle
if isMatrix
    A = get_or_apply_system_matrix(N,s,p,isDisp);
else
    A = @(U,TRANSP_FLAG) get_or_apply_system_matrix(...
        N,s,p,isDisp,U,TRANSP_FLAG);
end

% Create the phantom.
if nargout > 1
    x = phantomgallery('tectonic',N);
    x = x(:);
    if isMatrix
        b = A*x;
    else
        b = A(x,'notransp');
    end
end


function A = get_or_apply_system_matrix(N,s,p,isDisp,u,transp_flag)

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

% Prepare for illustration.
if isDisp
    AA = phantomgallery('tectonic',N); % rand(N);
    figure
end

% Deduce whether to set up matrix or apply to input, depending on whether
% input u is given.
isMatrix = (nargin < 5);

if isMatrix
    
    % Initialize the index for storing the results and the number of angles.
    idxend = 0;
    rows = zeros(2*N*s*p,1);
    cols = rows;
    vals = rows;
    
    II = 1:s;
    JJ = 1:p;
else
    % If transp_flag is 'size', only return size of operator, otherwise set
    % proper size of output.
    switch transp_flag
        case 'size'
            A = [p*s,N^2];
            return
        case 'notransp' % Forward project.
            if length(u) ~= N^2, error('Incorrect length of input vector'), end
            A = zeros(p*s,1);
        case 'transp' % Backproject
            if length(u) ~= p*s, error('Incorrect length of input vector'), end
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
        II = 1:s;
        JJ = 1:p;
    end
end

% Loop over all the sources.
for i = II
        
    % Illustation of the domain
    if isDisp

        clf
        pause(isDisp)
        imagesc((-N2+.5):(N2-0.5),(-N2+.5):(N2-0.5),flipud(AA))
        colormap gray
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
    for j = JJ
        
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
                                             
            set(gca,'Xtick',[],'Ytick',[])
            pause(isDisp)
        end
        
        % Collect the intersection coordinates and times.
        t = [tx; ty];
        xxy = [x; xy];
        yxy = [yx; y];
        
        % Sort the times.
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
        aval = sqrt(diff(xxy).^2 + diff(yxy).^2);
        %numvals = numel(aval);
        col = [];
        
        % Store the values inside the box.
        if numel(aval) > 0
                           
            % Calculates the midpoints of the line within the cells.
            xm = 0.5*(xxy(1:end-1) + xxy(2:end)) + N2;
            ym = 0.5*(yxy(1:end-1) + yxy(2:end)) + N2;
            
            % Translate the midpoint coordinates to index.
            col = floor(xm)*N + (N - floor(ym));

        end
        
        if ~isempty(col)
            if isMatrix
                % Create the indices to store the values to a vector for
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
    A = sparse(rows,cols,vals,s*p,N^2);
end
