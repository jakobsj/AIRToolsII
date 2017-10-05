function [A,b,x,angles,numCircles] = ...
    sphericaltomo(N,angles,numCircles,isDisp,isMatrix)
%SPHERICALTOMO  Creates a 2D spherical Radon tomography test problem
%
%   [A,b,x,angles,numCircles] = sphericaltomo(N)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles,numCircles)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles,numCircles,isDisp)
%   [A,b,x,angles,numCircles] = sphericaltomo(N,angles,numCircles,...
%                                             isDisp,isMatrix)
%
% This function genetates a tomography test problem based on the spherical
% Radon tranform where data consists of integrals along circles.  This type
% of problem arises, e.g., in photoacoustic imaging.
%
% The image domain is a square centered at the origin.  The centers for the
% integration circles are placed on a circle just outside the image domain.
% For each circle center we integrate along a number of concentric circles
% with equidistant radii, using the periodic trapezoidal rule.
%
% Input:
%   N           Scalar denoting the number of pixels in each dimesion, such
%               that the image domain consists of N^2 cells.
%   angles      Vector containing the angles to the circle centers in
%               degrees. Default: angles = 0:2:358.
%   numCircles  Number of concentric integration circles for each center.
%               Default: numCircles = round(sqrt(2)*N).
%   isDisp      If isDisp is non-zero it specifies the time in seconds
%               to pause in the display of the rays. If zero (the default),
%               no display is shown.
%   isMatrix    If non-zero, a sparse matrix is returned in A (default).
%               If zero, instead a function handle is returned.
%
% Output:
%   A           If input isMatrix is 1 (default): coefficient matrix with
%               N^2 columns and length(angles)*p rows.
%               If input isMatrix is 0: A function handle representing a
%               matrix-free version of A in which the forward and backward
%               operations are computed as A(x,'notransp') and A(y,'transp'),
%               respectively, for column vectors x and y of appropriate size.
%               The size of A can be retrieved using A([],'size'). The matrix
%               is never formed explicitly, thus saving memory.
%   b           Vector containing the rhs of the test problem.
%   x           Vector containing the exact solution, with elements
%               between 0 and 1.
%   angles      Similar to input.
%   numCircles  Similar to input.
%
% See also: paralleltomo, fancurvedtomo, fanlineartomo, seismictomo,
%           seismicwavetomo.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.
% Based on Matlab code written by Juergen Frikel, OTH Regensburg.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Default value of the angles to the circle centers.
if nargin < 2 || isempty(angles)
    angles = 0:2:358;
end

% Make sure angles is double precison to prevent round-off issues caused by
% single input.
angles = double(angles);

% Default value of the number of circles.
if nargin < 3 || isempty(numCircles)
    numCircles = round(sqrt(2)*N);
end

% Default illustration: False
if nargin < 4 || isempty(isDisp)
    isDisp = 0;
end

% Default matrix or matrix-free function handle? Matrix
if nargin < 5 || isempty(isMatrix)
    isMatrix = 1;
end

% Construct either matrix or function handle
if isMatrix
    A = get_or_apply_system_matrix(N,angles,numCircles,isDisp);
else
    A = @(U,TRANSP_FLAG) get_or_apply_system_matrix(N,angles,numCircles,...
        isDisp,U,TRANSP_FLAG);
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


function A = ...
    get_or_apply_system_matrix(N,angles,numCircles,isDisp,u,transp_flag)

% Define the number of angles.
nA = length(angles);

% Radii for the circles.
radii  = linspace(0,2,numCircles+1);
radii  = radii(2:end);

% Image coordinates.
centerImg = ceil(N/2);

% Determine the quarature parameters.
dx = sqrt(2)/N;
nPhi = ceil((4*pi/dx)*radii);
dPhi = 2*pi./nPhi;

% Prepare for illustration
if isDisp
    AA = phantomgallery('smooth',N);
    figure
end

% Deduce whether to set up matrix or apply to input, depending on whether
% input u is given.
isMatrix = (nargin < 5);

if isMatrix
    
    % Initialize vectors that contains the row numbers, the column numbers
    % and the values for creating the matrix A effiecently.
    rows = zeros(2*N*nA*numCircles,1);
    cols = rows;
    vals = rows;
    idxend = 0;
    
    II = 1:nA;
    JJ = 1:numCircles;
else
    % If transp_flag is 'size', only return size of operator, otherwise set
    % proper size of output.
    switch transp_flag
        case 'size'
            A = [numCircles*nA,N^2];
            return
        case 'notransp' % Forward project.
            if length(u) ~= N^2
                error('Incorrect length of input vector')
            end
            A = zeros(numCircles*nA,1);
        case 'transp' % Backproject
            if length(u) ~= numCircles*nA
                error('Incorrect length of input vector')
            end
            A = zeros(N^2,1);
    end
    
    % If u is a Cartesian unit vector then we only want to compute a single
    % row of A; otherwise we want to multiply with A or A'.
    if strcmpi(transp_flag,'transp') && nnz(u) == 1 && sum(u) == 1
        % Want to compute a single row of A, stored as a column vector.
        ell = find(u==1);
        JJ = mod(ell,numCircles);  if JJ==0, JJ = numCircles; end
        II = (ell-JJ)/numCircles + 1;
    else
        % Want to multiply with A or A'.
        II = 1:nA;
        JJ = 1:numCircles;
    end
end

% Loop over angles.
for m = II
    
    % Illustration of the domain.
    if isDisp
        clf
        pause(isDisp)
        imagesc(0.5:N-0.5,0.5:N-0.5,flipud(AA))
        colormap gray
        axis xy
        hold on
        axis equal
        axis([-0.3*N 1.3*N -0.3*N 1.3*N])
    end
    
    % Angular position of source.
    xix = cosd(angles(m));
    xiy = sind(angles(m));
    
    % Loop over the circles.
    for n = JJ
        
        % (x,y) coordinates of circle.
        k = (0:nPhi(n))*dPhi(n);
        xx = (xix + radii(n)*cos(k))/dx + centerImg;
        yy = (xiy + radii(n)*sin(k))/dx + centerImg;
        
        % Illustration of the circle.
        if isDisp
            
            plot(xx,yy,'-','color',[220 0 0]/255,'linewidth',1.5)
            
            set(gca,'Xtick',[],'Ytick',[])
            pause(isDisp)
        end
        
        % Round to get pixel index.
        col = round( xx );
        row = round( yy );
        
        % Discard if outside box domain.
        IInew = find(col>0 & col<=N & row>0 & row<=N);
        row = row(IInew);
        col = col(IInew);
        J = (N-row+1) + (col-1)*N;
        
        % Convert to linear index and bin
        Ju = unique(J);
        w = histc(J,Ju);
        
        % Determine rows, columns and weights.
        i = (m-1)*numCircles + n;
        ii = repmat(i,length(Ju),1);
        jj = Ju(:);
        aa = (2*pi*radii(n)/nPhi(n))*w(:);
        
        % Store the values, if any.
        if ~isempty(jj)
            if isMatrix
                % Create the indices to store the values to vector for
                % later creation of A matrix.
                idxstart = idxend + 1;
                idxend = idxstart + numel(jj) - 1;
                idx = idxstart:idxend;
                
                % Store row numbers, column numbers and values.
                rows(idx) = ii;
                cols(idx) = jj;
                vals(idx) = aa;
            else
                % If any nonzero elements, apply forward or back operator
                
                if strcmp(transp_flag,'notransp')
                    % Insert inner product with u into w.
                    A( i ) = aa'*u(jj);
                else  % Adjoint operator.
                    A(jj) = A(jj) + u( i )*aa;
                end
            end
        end
    end
end

if isMatrix
    % Truncate excess zeros.
    rows = rows(1:idxend);
    cols = cols(1:idxend);
    vals = vals(1:idxend);
    
    % Create sparse matrix A from the stored values.
    A = sparse(rows,cols,vals,numCircles*nA,N^2);
end