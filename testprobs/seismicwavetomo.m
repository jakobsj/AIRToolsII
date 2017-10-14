function [A,b,x,s,p,omega] = seismicwavetomo(N,s,p,omega,isDisp,isMatrix)
%SEISMICWAVETOMO  Seismic tomography problem without the ray assumption
%
%   [A,b,x,s,p] = seismicwavetomo(N)
%   [A,b,x,s,p] = seismicwavetomo(N,s)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p,omega)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p,omega,isDisp)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p,omega,isDisp,isMatrix)
%
% This function creates a 2D seismic travel-time tomography test problem
% with an N-times-N pixel domain, using s sources located on the right
% boundary and p receivers (seismographs) scattered along the left and top
% boundary. Waves are transmitted from each source to each receiver, and
% there is no high-frequency approximation, i.e., no assumption about rays.
% The wave is assumed to travel within the first Fresnel zone.
%
% In the limit where the wave frequency is infinite, the waves can be
% described by rays as done in the function seismictomo. But note that
% A and b do not converge to those of seismictomo as omega -> inf.
%
% Inupt: 
%   N         Scalar denoting the number of discretization intervals in 
%             each dimesion, such that the domain consists of N^2 cells.
%   s         The number of sources in the right side of the domain
%             (default s = N.
%   p         The number of receivers (seismographs) equally spaced on the
%             surface and on the left side of the domain. Default p = 2*N.
%   omega     Dominant frequency of the propagating wave. Default = 10.
%   isDisp    If isDisp is nonzero it specifies the time in seconds to pause
%             in the display of the rays. If zero (the default), then no
%             display is shown.
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
%   omega     The frequency.
%
% See also: paralleltomo, fancurvedtomo, fanlineartomo, seismictomo,
%           sphericaltomo.

% Reference: J. M. Jensen, B. H. Jacobsen, and J. Christensen-Dalsgaard,
% Sensitivity kernels for time-distance inversion, Solar Physics, 192
% (2000), pp. 231-239.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.
% With contribution from Mikkel Brettschneider.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Default number of sources.
if nargin < 2 || isempty(s)
    s = N;
end

% Default number of receivers (seismographs).
if nargin < 3 || isempty(p)
    p = 2*N;
end

% Default frequency.
if nargin < 4 || isempty(omega)
    omega = 10;
end

if nargin < 5 || isempty(isDisp)
    isDisp = 0;
end

if nargin < 6 || isempty(isMatrix)
    isMatrix = 1;
end

% Construct either matrix or function handle
if isMatrix
    A = get_or_apply_system_matrix(N,s,p,omega,isDisp);
else
    A = @(U,TRANSP_FLAG) get_or_apply_system_matrix(...
        N,s,p,omega,isDisp,U,TRANSP_FLAG);
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



function A = get_or_apply_system_matrix(N,s,p,omega,isDisp,u,transp_flag)

% Threshold for sparsity.
tol = 1e-6;

N2 = N/2;

% Determine the positions of the sources.
Ns = (N/s)/2;
x0 = N2*ones(s,1);                 
y0 = linspace(-N2+Ns,N2-Ns,s)';

% The positions of the seismographs.
p1 = ceil(p/2); 	
p2 = floor(p/2); 
Np1 = (N/p1)/2;                 
Np2 = (N/p2)/2;
xp = [-N2*ones(p2,1); linspace(-N2+Np1,N2-Np1,p1)'];
yp = [linspace(-N2+Np2,N2-Np2,p2)'; N2*ones(p1,1)];

% The intersection lines.
xrange = (-N2+.5:N2-.5)';
yrange = (N2-.5:-1:-N2+.5)';
[xx,yy] = meshgrid(xrange,yrange);

% Prepare for illustration.
if isDisp
    figure
end

% Deduce whether to set up matrix or apply to input from whether input u is
% given.
isMatrix = (nargin < 6);

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
        pause(isDisp)
        clf;
        hold on
        axis xy equal
        axis([-N2 N2 -N2 N2])
        plot(x0,y0,'.','color','r','linewidth',2,'markersize',30)
        plot(xp,yp,'s','MarkerEdgeColor','m',...
                'MarkerFaceColor','m','MarkerSize',10,'linewidth',2)
    end
        
    % Loop over the seismographs.
    for j = JJ
        
        % Kreate sensitivity kernels
        R = [xp(j) yp(j)];
        S = [x0(i) y0(i)];
        sens = simple_fat_kernel(R,S,omega/N,xx,yy);
        sens(sens<tol) = 0;

        [ym,xm,aval] = find(sens);

        % Illustration of the sensitivity kernels.
        if isDisp
            imagesc(xrange,yrange,sens);
            set(gca,'Xtick',[],'Ytick',[])
            pause(isDisp)
        end

        %numvals = numel(aval);
        col = [];
        
        % Store the values inside the domain.
        if numel(aval) > 0
                                       
            % Translate the domain index to matrix index.
            col = (xm-1)*N + ym;

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



function S = simple_fat_kernel(R,S,omega,xx,yy)
% Computation of Fresnel zones that describe the sensitivity of the solution
% to the seismic wave.

alpha = 10;  % Chosen by inspection.

tSX = sqrt((xx-S(1)).^2 + (yy-S(2)).^2);
tRX = sqrt((xx-R(1)).^2 + (yy-R(2)).^2);
distSR = sqrt(sum((S-R).^2));
delta_t = tSX+tRX - distSR;

% Compute the sensitivity kernel.
S = cos(2*pi*delta_t.*omega).*exp(- (alpha*delta_t.*omega).^2 );

% Normalize S.
S = distSR.*S./sum(S(:));
