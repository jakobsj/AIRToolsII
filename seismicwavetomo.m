function [A,b,x,s,p,omega] = seismicwavetomo(N,s,p,omega,isDisp)
%SEISMICWAVETOMO  Seismic tomography problem without the ray assumption
%
%   [A,b,x,s,p] = seismicwavetomo(N)
%   [A,b,x,s,p] = seismicwavetomo(N,s)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p,omega)
%   [A,b,x,s,p] = seismicwavetomo(N,s,p,omega,isDisp)
%
% This function creates a 2D seismic tomography test problem with an
% N-times-N domain, using s sources located on the right boundary and p
% receivers (seismographs) scattered along the left and top boundary.
% Waves are transmitted from each source to each receiver, and there is
% no high-frequency approximation, i.e., no assumption about rays. The
% wave is assumed to travel within the first Fresnel zone.
%
% In the limit where the wave frequency is infinite, the waves can be
% described by rays as done in the function seismictomo. But note that
% A and b do not converge to those of seismictomo as omega -> inf.
%
% Inupt: 
%   N        Scalar denoting the number of discretization intervals in 
%            each dimesion, such that the domain consists of N^2 cells.
%   s        The number of sources in the right side of the domain
%            (default s = N.
%   p        The number of receivers (seismographs) equally spaced on the
%            surface and on the left side of the domain (default p = 2*N).
%   omega    Dominant frequency of the propagating wave (default = 10).
%   isDisp   If isDisp is nonzero it specifies the time in seconds to pause
%            in the display of the rays. If zero (the default), then no
%            display is shown.
%
% Output:
%   A        Coefficient matrix with N^2 columns and s*p rows.
%   b        Vector containing the rhs of the test problem.
%   x        Vector containing the exact solution, with elements
%            between 0 and 1.
%   s        The number of sources.
%   p        The number of receivers (seismographs).
%   omega    The frequency.
%
% See also: seismictomo.

% Mikkel Brettschneider and Per Chr. Hansen, December 15, 2012.
% Based on seismictomo from AIR Tools.

% Reference: J. M. Jensen, B. H. Jacobsen, and J. Christensen-Dalsgaard,
% Sensitivity kernels for time-distance inversion, Solar Physics, 192
% (2000), pp. 231-239.
        
% Threshold for sparsity.
tol = 1e-6;

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

% Initialize the arrays for creating the sparse matrix.
idxend = 0;
rows = zeros(N^2*s*p,1);
cols = rows;
vals = rows;

% Prepare for illustration.
if isDisp
    figure
end

% Loop over all the sources.
for i = 1:s
        
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
    for j = 1:p
        
        % Kreate sensitivity kernels
        R = [xp(j) yp(j)];
        S = [x0(i) y0(i)];
        sens = simple_fat_kernel(R,S,omega/N,xx,yy);
        sens(sens<tol) = 0;

        [ym,xm,d] = find(sens);

        % Illustration of the sensitivity kernels.
        if isDisp
            imagesc(xrange,yrange,sens);
            set(gca,'Xticklabel',{},'Yticklabel',{})
            pause(isDisp)
        end

        numvals = numel(d);
        
        % Store the values inside the domain.
        if numvals > 0
                                       
            % Translate the domain index to matrix index.
            col = (xm-1)*N + ym;
            
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

% Right-hand side and solution, if needed.
if nargout > 1
    x = tectonic(N); % Creates the phantom.
    b = A*x;
end

end

% Subfunctions start here --------------------------------------------

function x = tectonic(N)
% Creates a tectonic phantom of size N-by-N.

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

end