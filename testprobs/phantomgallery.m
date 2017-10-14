function im = phantomgallery(name,N,P1,P2,P3)
%PHANTOMGALLERY  A collection of 2D phantoms for use in test problems
%
%   im = phantomgallery(name,N)
%   im = phantomgallery(name,N,P1)
%   im = phantomgallery(name,N,P1,P2)
%   im = phantomgallery(name,N,P1,P2,P3)
%
% The phantom im is N-by-N with pixel values between 0 and 1, and the
% following phantom types are available:
%
% shepplogan: the Shepp-Logan phantom
%   im = phantomgallery('shepplogan',N)
%
% smooth: a smooth image
%   im = phantomgallery('smooth',N,P1)
%   P1 = 1, 2, 3 or 4 defines four different smooth functions (default = 4)
% The image is constructed by adding four different Gaussian functions.
%
% binary: a random image with binary pixel values arranged in domains
%   im = phantomgallery('binary',N,P1)
%   P1 = seed for random number generator
% The image is dominated by horizontal structures.
%
% threephases: a random image with pixel values 0, 0.5, 1 arranged in domains
%   im = phantomgallery('threephases',N,P1,P2)
%   P1 = controls the number of domains (default = 100)
%   P2 = seed for random number generator
% The image is a model of a three-phase object.
%
% threephasessmooth: similar to threephases, but the domains have smoothly
%                    varying pixel values and there is a smooth background
%   im = phantomgallery('threephasessmooth',N,P1,P2,P3)
%   P1 = controls the number of domains (default = 100)
%   P2 = controls the intensity variation within each domain (default = 1.5)
%   P3 = seed for random number generator
%
% fourphases: a random image similar to 'binary' but with three phases
%             separated by (thin) structures that form the fourth phase
%   im = phantomgallery('fourphases',N,P1)
%   P1 = seed for random number generator
%
% grains: a random image with Voronoi cells
%   im = phantomgallery('grains',N,P1,P2)
%   P1 = number of cells in the image (default = 3*sqrt(N))
%   P2 = seed for random number generator
% The image is a model of grains with different pixel intensities.
%
% ppower: a random image with patterns of nonzero pixels
%   im = phantomgallery('ppower',N,P1,P2,P3)
%   P1 = the ratio of nonzero pixels, between 0 and 1 (default = 0.3)
%   P2 = the smoothness of the image, greater than 0 (default = 2)
%   P3 = seed for random number generator
% The larger the P2 the larger the domains of nonzero pixels.
%
% tectonic: a test image for the seismic tomography test problems.
%   im = phantomgallery('tectonic',N)
%
% To use these images in connection with the test problems use the commands:
%   im = phantomgallery(name,N,...);
%   x = im(:);
%   A = matrix generated, e.g., by paralleltomo;
%   b = A*x;

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.
% With contibutions by Mikhail Romanov, Technical University of Denmark
% and Knud Cordua, University of Copenhagen.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

if nargin < 2, error('Not enought input arguments'), end

switch name
    
    case 'shepplogan'
        im = myphantom(N);
    
    case 'smooth'
        if nargin == 2
            im = smooth(N);
        else
            im = smooth(N,P1);
        end
        
    case 'binary'
        if nargin == 2
            im = binary(N);
        else
            im = binary(N,P1);
        end
        
    case 'threephases'
        if nargin == 2
            im = threephases(N);
        elseif nargin == 3
            im = threephases(N,P1);
        else
            im = threephases(N,P1,P2);
        end
        
    case 'threephasessmooth'
        if nargin ==2
            im = threephasessmooth(N);
        elseif nargin == 3
            im = threephasessmooth(N,P1);
        elseif nargin == 4
            im = threephasessmooth(N,P1,P2);
        else
            im = threephasessmooth(N,P1,P2,P3);
        end
        
    case 'fourphases'
        if nargin == 2
            im = fourphases(N);
        else
            im = fourphases(N,P1);
        end
        
    case 'grains'
        if nargin == 2
            im = grains(N);
        elseif nargin == 3
            im = grains(N,P1);
        else
            im = grains(N,P1,P2);
        end
        
    case 'ppower'
        if nargin ==2
            im = ppower(N);
        elseif nargin == 3
            im = ppower(N,P1);
        elseif nargin == 4
            im = ppower(N,P1,P2);
        else
            im = ppower(N,P1,P2,P3);
        end
        
    case 'tectonic'
        im = tectonic(N);
        
    otherwise
        error('Illegal phantom name')
        
end

% Subfunctions here ----------------------------------------------------

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

% Ensure nonnegative elements.
X(X<0) = 0;

function im = smooth(N,p)
%SMOOTH Creates a 2D test image of a smooth function

% Per Christian Hansen, May 8, 2012, DTU Compute.

if nargin==1, p = 4; end

% Generate the image.
[I,J] = meshgrid(1:N);
sigma = 0.25*N;
c = [0.6*N 0.6*N; 0.5*N 0.3*N; 0.2*N 0.7*N; 0.8*N 0.2*N];
a = [1 0.5 0.7 0.9];
im = zeros(N,N);
for i=1:p
    im = im + a(i)*exp( - (I-c(i,1)).^2/(1.2*sigma)^2 - (J-c(i,2)).^2/sigma^2);
end
im = im/max(im(:));

% -----------------------------------------------------------------------

function im = binary(N,seed)
%BINARY Creates a 2D binary test image

% Per Christian Hansen, June 9, 2013, DTU Compute.  This function uses
% software kindly provided by Knud Cordua, Univ. of Copenhagen.

if nargin==2, rng(seed), end

% Prepare to generate random image.
TrainingImage = getTI;

C = ones(3);  % Needed in the two functions below.
nCat = 1;     % Ditto.
H = CompHistTrain(TrainingImage,C,nCat) + ...
    CompHistTrain(fliplr(TrainingImage),C,nCat) + eps;

% Generate a random image.
im = sequential_simulation_full_clique(H,nCat,N,N,C);

% -----------------------------------------------------------------------

function im = threephases(N,p,seed)
%THREEPHASES Creates a 2D test image with three different phases

% Per Christian Hansen, Sept. 30, 2014, DTU Compute.

if nargin==1 || isempty(p), p = 100; end
if nargin==3, rng(seed), end

% Generate first image.
[I,J] = meshgrid(1:N);
sigma1 = 0.025*N;
c1 = rand(p,2)*N;
im1 = zeros(N,N);
for i=1:p
    im1 = im1 + exp(-abs(I-c1(i,1)).^3/(2.5*sigma1)^3 ...
                    -abs(J-c1(i,2)).^3/sigma1^3);
end
t1 = 0.35;
im1(im1 < t1) = 0;
im1(im1 >= t1) = 2;

% Generate second image.
sigma2 = 0.025*N;
c2 = rand(p,2)*N;
im2 = zeros(N,N);
for i=1:p
    im2 = im2 + exp(-(I-c2(i,1)).^2/(2*sigma2)^2-(J-c2(i,2)).^2/sigma2^2);
end
t2 = 0.55;
im2(im2 < t2) = 0;
im2(im2 >= t2) = 1;

% Combine the two images.
im = im1 + im2;
im(im == 3) = 1;
im = im/max(im(:));

% -----------------------------------------------------------------------

function im = threephasessmooth(N,p,v,seed)
%THREEPHASESSMOOTH Variant of threephases with smoothly varying intensities.

% Per Christian Hansen and Mikhail Romanov, July 6, 2015, DTU Compute.

if nargin==1 || isempty(p), p = 100; v = 1.8; end
if nargin==2, v = 1.8; end
if nargin==4, rng(seed), end

% Generate first image.
[I,J] = meshgrid(1:N);
sigma1 = 0.025*N;
c1 = rand(p,2)*N;
im1 = zeros(N,N);
for i=1:p
    im1 = im1 + exp(-abs(I-c1(i,1)).^3/(2.5*sigma1)^3 ...
                    -abs(J-c1(i,2)).^3/sigma1^3);
end
t1 = 0.35;
im1(im1 < t1) = 0;
I1 = find(im1 >= t1);
im1(I1) = (im1(I1) - min(im1(I1)))/max(im1(:))*v + 0.8;

% Generate second image.
sigma2 = 0.025*N;
c2 = rand(p,2)*N;
im2 = zeros(N,N);
for i=1:p
    im2 = im2 + exp(-(I-c2(i,1)).^2/(2*sigma2)^2-(J-c2(i,2)).^2/sigma2^2);
end
t2 = 0.55;
im2(im2 < t2) = 0;
I2 = find(im2 >= t2);
im2(I2) = (im2(I2) - min(im2(I2)))/max(im2(:))*v + 0.3;

% Combine the two images onto a smooth background..
im = (v/3)*ppower(N,1,2.5);
im(im1 > 0) = im1(im1 > 0);
im(im2 > 0) = im2(im2 > 0);
im = im/max(im(:));

% -----------------------------------------------------------------------

function im = fourphases(N,seed)
%FOURPHASES Creates a 2D test image with three different phases

% Per Christian Hansen and Mikhail Romanov, March. 1, 2015, DTU Compute.

if nargin==2, rng(seed), end

x = 1 - phantomgallery('binary',N);
x = bwlabel(x);
im = (mod(x+1,4) == 0)*0.33 + (mod(x+2,4) == 0)*0.66 + (mod(x+3,4) == 0);

% -----------------------------------------------------------------------

function im = grains(N,numGrains,seed)
%GRAINS Creates a test image of Voronoi cells

% Jakob Sauer Jorgensen, October 9, 2012, DTU Compute.

if nargin==1 || isempty(numGrains), numGrains = round(3*sqrt(N)); end
if nargin==3, rng(seed), end

% Prepare to create a bigger image, to avoid strange boundary grains.
dN = round(N/10);
Nbig = N + 2*dN;
total_dim = Nbig^2;

% Random pixels whose coordinates (xG,yG,zG) are the "centre" of the grains.
xG = ceil(Nbig*rand(numGrains,1));
yG = ceil(Nbig*rand(numGrains,1));

% Set up voxel coordinates for distance computation.
[X,Y] = meshgrid(1:Nbig);
X     = X(:);
Y     = Y(:);

% For centre pixel k [xG(k),yG(k),zG(k)] compute the distance to all the 
% voxels in the box and store the distances in column k.
distArray = zeros(total_dim,numGrains);
for k = 1:numGrains
    distArray(:,k) = (X-xG(k)).^2 + (Y-yG(k)).^2;
end

% Determine to which grain each of the voxels belong. This is found as the
% centre with minimal distance to the given voxel.
[~,minIdx] = min(distArray,[],2);

% Reshape to 2D, subtract 1 to have 0 as minimal value, extract the
% middle part of the image, and scale to have 1 as maximum value.
im = reshape(minIdx,repmat(Nbig,1,2)) - 1;
im = im(dN+(1:N),dN+(1:N));
im = im/max(im(:));

% -----------------------------------------------------------------------

function F = ppower(N,relnz,p,seed)
%PPOWER Creates a 2D test image with patterns of nonzero pixels

% Per Christian Hansen and Jakob Sauer Jorgensen, July 6, 2015.

% Reference: J. S. Jorgensen, E. Y. Sidky, P. C. Hansen, and X. Pan,
% Empirical average-case relation between undersampling and sparsity in
% X-ray CT, submitted.

if nargin<2 || isempty(relnz), relnz = 0.3; end
if nargin<3 || isempty(p), p = 2; end
if nargin==4, rng(seed), end

if N/2 == round(N/2), Nodd = false; else Nodd = true; N = N+1; end

P = randn(N,N);
[I,J] = meshgrid(1:N);
U = ( ( (2*I-1)/N - 1).^2 + ( (2*J-1)/N - 1).^2 ).^(-p/2);
F = U.*exp(2*pi*sqrt(-1)*P);
F = abs(ifft2(F));
f = sort(F(:),'descend');
k = round(relnz*N^2);
F( F < f(k) ) = 0;
F = F/f(1);

if Nodd, F = F(1:end-1,1:end-1); end

% -----------------------------------------------------------------------

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

% ------ Below here: subfunctions for 'binary' --------------------------

function H = CompHistTrain(Z, N, sV)
%
% INPUTS
%  Z   binary valued training image (for now... )
%  N   binary neighborhood mask
%  sV  the image must have the categories 0,1,...,sV
%
% OUTPUT
%  H   image

nc = 1; mc = 1; pc = 1;

% Size of the training image and the neighbourhood mask.
[nZ, mZ, pZ] = size(Z);
[nN, mN, pN] = size(N);
sN = sum(N(:));

% Base of the conversion from neighborhood to ID
base = ((sV+1)*ones(1,sN)).^(0:sN-1);

% Indexes of the voxels
index = reshape(1:nZ*mZ*pZ, [nZ, mZ, pZ]);

% Histogram
H = zeros(1,(sV+1)^sN);

% For all inner voxels
for k = pc : pZ-(pN-pc)
   for j = mc : mZ-(mN-mc)
      for i = nc : nZ-(nN-nc)
         
         % The indexes of the rectangular neighborhood
         dim1 = i-nc+1 : i-nc+nN; 
         dim2 = j-mc+1 : j-mc+mN; 
         dim3 = k-pc+1 : k-pc+pN; 

         % Neighbors ID of the voxel
         ijk = index(dim1,dim2,dim3);
                  
         % Voxel values of Nk
         zijk = Z(ijk);
         nijk = zijk(N==1);
       
         % Identify neighborhood id:
         hijk = n2id(nijk(:), base);
         
         
         %Modified by ksc:
         H(1, hijk) = H(1, hijk) + 1;
         
      end
   end
end
 
function sim = sequential_simulation_full_clique(H,sV,ny,nx,C)
% 
% Input parameters:
%   H:  Histogram of patterns
%   Sv: Number of categories (Sv=1 for two categories).
%   ny: Number of vertical nodes in the simulated field
%   nx: Number of horizontal nodes in the simulated field
%   C:  The clique (i.e. template for pattern scanning)
%   N_real: Number of realizations to be produced.
%
% Output parameter:
% sim:  An array with size = [ny nx] that contains the realization.

H_full_clique = H;

% Size of clique;
[NyC,NxC] = size(C);

% Initialize arrays:
sim  = ones(ny,nx).*NaN;
prob = zeros(ny,nx);

% Simulate the first unconditional pattern:
norm = sum(H_full_clique);
cdf  = cumsum(H_full_clique)/norm;
r    = rand;
id   = find(cdf>r,1);
prob(1,1) = H_full_clique(id)/norm;
sN = NyC*NxC;
n  = id2n(id,sV,sN); % Find the pattern associated with the id
sim(1:NyC,1:NxC) = reshape(n,NyC,NxC);

    % Simulate the first row of residuals/seperators:
    i=1:NyC;
    for j=NxC+1:nx
        
        % Position of the full clique to be considered:
        marginal=sim(i,j-NxC+1:j);
        marginal(isnan(marginal))=-1;
        % Find the position at which the known values are located:
        ids=marg2id(marginal(:));
        pdf=H_full_clique(ids)/sum(H_full_clique(ids));
        cdf=cumsum(pdf)/sum(pdf);
        r=rand;
        pos=find(cdf>r,1);
        id=ids(pos);
        prob(1,j)=H_full_clique(id)/sum(H_full_clique);
        n=id2n(id,sV,sN); % Find the pattern associated with the id
        n=reshape(n,NyC,NxC);
        sim(i,j)=n(:,end);
        
    end
    
    for i=NyC+1:ny
        for j=NxC:nx
            
            if j==NxC % Left side of the simulation area
                % Position of the full clique to be considered:
                marginal=sim(i-NyC+1:i,j-NxC+1:j);
                marginal(isnan(marginal))=-1;
                ids=marg2id(marginal(:));
                
                if sum(H_full_clique(ids))==0,
                    error('Zero probabiltiy pattern event! Please add a homogeneous distribution'),
                end
                
                pdf=H_full_clique(ids)/sum(H_full_clique(ids));
                cdf=cumsum(pdf)/sum(pdf);
                r=rand;
                pos=find(cdf>r,1);
                id=ids(pos);
                prob(i,j)=H_full_clique(id)/sum(H_full_clique);
                n=id2n(id,sV,sN); % Find the pattern associated with the id
                n=reshape(n,NyC,NxC);
                sim(i,1:j)=n(end,:);
                
            else % Rest of the nodes to be simulated:
                
                % Position of the full clique to be considered:
                marginal=sim(i-NyC+1:i,j-NxC+1:j);
                marginal(isnan(marginal))=-1;
                ids=marg2id(marginal(:));
                
                if sum(H_full_clique(ids))==0,
                    error('Zero probabiltiy pattern event! Please add a homogeneous distribution'),
                end
                
                pdf=H_full_clique(ids)/sum(H_full_clique(ids));
                cdf=cumsum(pdf)/sum(pdf);
                r=rand;
                pos=find(cdf>r,1);
                id=ids(pos);
                prob(i,j)=H_full_clique(id)/sum(H_full_clique);
                n=id2n(id,sV,sN);
                n=reshape(n,NyC,NxC);
                sim(i,j)=n(end,end);
                
            end
        end
    end
    
function n = id2n(id, sV, sN)
%
% INPUT
%     id    scalar, integer ID number of the neighborhood, id > 0
%     sV    scalar, number of different values the voxels can take, sV = 1
%           for binary images
%     sN    scalar, number of neighbors of a voxel (not on the boundary)
%
% OUTPUT
%     n     vector, sN x 1, with the values of the voxels in the neighborhood

n = zeros(sN,1);
id = id - 1;

for i = 1:sN
   idi = floor(id / (sV+1));
   n(i) = id - (sV+1)*idi;
   id = idi;
end

function id = n2id(n, base)
% n2id: assigns ID number to a neighborhood
%
% INPUT
%     n     vector, sN x 1, with the values of the voxels in the
%           neighborhood
%     base  vector, 1 x sN, base = (sV*ones(1,sN)).^(0:sN-1), see id2n for
%           definitions of sV and sN
%
% OUTPUT
%     id    scalar, integer ID number of the neighborhood, id > 0.

id = base * n + 1;

function id=marg2id(marginal)

% The id outcomes of this function should be used to calculate marginals by 
% summing the probabilities related to these ids. 

mm = zeros(length(marginal),1);
mm(marginal>-1) = 1;
marg_dim = sum(mm);

sV = 2; %(only tested for binary images, sV=2)
sN = length(marginal);

index = 1:sV^sN;
Index = zeros(sV^sN,marg_dim); 

c = 0;
for i=1:sN
    if marginal(i) >= 0
        c = c+1;
        % All indeces:
        index_tmp = reshape(index,sV^(i-1),sV^sN/sV^(i-1));

        % Indeces used for the particular numbers:  
        index_use = index_tmp(:,marginal(i)+1:2:sV^sN/sV^(i-1));
        Index(index_use(:),c) = 1;
    end
end

id = find(sum(Index,2)==marg_dim);

function A = getTI
%
% Create the binary training image.

IJ = [
9 1
10 1
22 1
23 1
24 1
46 1
47 1
48 1
59 1
60 1
61 1
76 1
77 1
78 1
79 1
9 2
10 2
22 2
23 2
24 2
25 2
45 2
46 2
47 2
48 2
59 2
60 2
61 2
75 2
76 2
77 2
78 2
8 3
9 3
10 3
22 3
23 3
24 3
25 3
26 3
44 3
45 3
46 3
47 3
58 3
59 3
60 3
74 3
75 3
76 3
77 3
8 4
9 4
10 4
21 4
22 4
23 4
25 4
26 4
27 4
44 4
45 4
58 4
59 4
60 4
73 4
74 4
75 4
76 4
8 5
9 5
10 5
21 5
22 5
23 5
26 5
27 5
28 5
43 5
44 5
45 5
57 5
58 5
59 5
72 5
73 5
74 5
75 5
8 6
9 6
10 6
20 6
21 6
22 6
28 6
29 6
30 6
42 6
43 6
44 6
55 6
56 6
57 6
58 6
59 6
71 6
72 6
73 6
8 7
9 7
10 7
20 7
21 7
29 7
30 7
31 7
42 7
43 7
53 7
54 7
55 7
56 7
57 7
58 7
59 7
70 7
71 7
72 7
73 7
8 8
9 8
10 8
20 8
21 8
29 8
30 8
31 8
42 8
43 8
53 8
54 8
55 8
58 8
59 8
69 8
70 8
71 8
72 8
8 9
9 9
10 9
19 9
20 9
21 9
30 9
31 9
42 9
43 9
52 9
53 9
54 9
58 9
59 9
69 9
70 9
71 9
8 10
9 10
10 10
19 10
20 10
21 10
30 10
31 10
32 10
41 10
42 10
43 10
51 10
52 10
53 10
58 10
59 10
60 10
68 10
69 10
70 10
8 11
9 11
10 11
19 11
20 11
21 11
30 11
31 11
32 11
41 11
42 11
43 11
50 11
51 11
52 11
58 11
59 11
60 11
68 11
69 11
70 11
9 12
10 12
11 12
19 12
20 12
21 12
30 12
31 12
32 12
41 12
42 12
43 12
49 12
50 12
51 12
59 12
60 12
68 12
69 12
70 12
10 13
11 13
12 13
19 13
20 13
21 13
30 13
31 13
32 13
41 13
42 13
43 13
49 13
50 13
51 13
59 13
60 13
68 13
69 13
70 13
10 14
11 14
12 14
13 14
19 14
20 14
21 14
30 14
31 14
32 14
41 14
42 14
43 14
49 14
50 14
51 14
59 14
60 14
61 14
68 14
69 14
70 14
10 15
11 15
12 15
13 15
14 15
15 15
19 15
20 15
21 15
30 15
31 15
32 15
40 15
41 15
42 15
43 15
49 15
50 15
51 15
60 15
61 15
62 15
68 15
69 15
70 15
11 16
12 16
13 16
14 16
15 16
16 16
17 16
20 16
21 16
30 16
31 16
32 16
40 16
41 16
42 16
43 16
49 16
50 16
51 16
60 16
61 16
62 16
69 16
70 16
71 16
11 17
12 17
15 17
16 17
17 17
18 17
19 17
20 17
21 17
30 17
31 17
32 17
39 17
40 17
41 17
42 17
43 17
44 17
49 17
50 17
51 17
61 17
62 17
63 17
70 17
71 17
72 17
11 18
12 18
17 18
18 18
19 18
20 18
21 18
22 18
31 18
32 18
33 18
38 18
39 18
40 18
42 18
43 18
44 18
49 18
50 18
51 18
61 18
62 18
63 18
71 18
72 18
11 19
12 19
19 19
20 19
21 19
22 19
31 19
32 19
33 19
38 19
39 19
43 19
44 19
50 19
51 19
62 19
63 19
64 19
71 19
72 19
73 19
11 20
12 20
13 20
20 20
21 20
22 20
31 20
32 20
33 20
37 20
38 20
39 20
43 20
44 20
50 20
51 20
52 20
62 20
63 20
64 20
72 20
73 20
12 21
13 21
20 21
21 21
22 21
31 21
32 21
33 21
37 21
38 21
43 21
44 21
45 21
50 21
51 21
52 21
63 21
64 21
65 21
72 21
73 21
74 21
12 22
13 22
14 22
21 22
22 22
23 22
24 22
32 22
33 22
36 22
37 22
38 22
44 22
45 22
50 22
51 22
52 22
63 22
64 22
65 22
73 22
74 22
13 23
14 23
22 23
23 23
24 23
32 23
33 23
36 23
37 23
44 23
45 23
50 23
51 23
52 23
63 23
64 23
65 23
73 23
74 23
75 23
13 24
14 24
15 24
23 24
24 24
25 24
32 24
33 24
34 24
35 24
36 24
37 24
44 24
45 24
46 24
51 24
52 24
63 24
64 24
73 24
74 24
75 24
14 25
15 25
16 25
24 25
25 25
32 25
33 25
34 25
35 25
36 25
44 25
45 25
46 25
51 25
52 25
62 25
63 25
64 25
73 25
74 25
75 25
76 25
1 26
14 26
15 26
16 26
25 26
26 26
32 26
33 26
34 26
35 26
44 26
45 26
46 26
51 26
52 26
53 26
62 26
63 26
64 26
73 26
74 26
75 26
76 26
1 27
15 27
16 27
17 27
25 27
26 27
32 27
33 27
44 27
45 27
46 27
51 27
52 27
53 27
61 27
62 27
63 27
72 27
73 27
74 27
75 27
76 27
77 27
1 28
2 28
15 28
16 28
17 28
25 28
26 28
32 28
33 28
34 28
44 28
45 28
46 28
51 28
52 28
53 28
61 28
62 28
63 28
72 28
73 28
74 28
76 28
77 28
1 29
2 29
15 29
16 29
25 29
26 29
31 29
32 29
33 29
34 29
35 29
45 29
46 29
52 29
53 29
60 29
61 29
62 29
72 29
73 29
76 29
77 29
78 29
1 30
2 30
15 30
16 30
25 30
26 30
30 30
31 30
32 30
33 30
34 30
35 30
36 30
45 30
46 30
52 30
53 30
54 30
60 30
61 30
71 30
72 30
73 30
76 30
77 30
78 30
1 31
2 31
3 31
14 31
15 31
16 31
24 31
25 31
26 31
30 31
31 31
32 31
34 31
35 31
36 31
44 31
45 31
46 31
47 31
52 31
53 31
54 31
59 31
60 31
61 31
71 31
72 31
73 31
77 31
78 31
2 32
3 32
13 32
14 32
15 32
24 32
25 32
30 32
31 32
35 32
36 32
43 32
44 32
45 32
46 32
47 32
48 32
53 32
54 32
59 32
60 32
71 32
72 32
77 32
78 32
2 33
3 33
4 33
13 33
14 33
15 33
24 33
25 33
29 33
30 33
31 33
35 33
36 33
42 33
43 33
44 33
46 33
47 33
48 33
53 33
54 33
59 33
60 33
71 33
72 33
77 33
78 33
79 33
3 34
4 34
13 34
14 34
24 34
25 34
29 34
30 34
35 34
36 34
37 34
41 34
42 34
43 34
46 34
47 34
48 34
53 34
54 34
55 34
59 34
60 34
70 34
71 34
72 34
77 34
78 34
79 34
3 35
4 35
12 35
13 35
14 35
24 35
25 35
29 35
30 35
35 35
36 35
37 35
41 35
42 35
43 35
46 35
47 35
48 35
54 35
55 35
56 35
58 35
59 35
60 35
70 35
71 35
72 35
77 35
78 35
79 35
3 36
4 36
5 36
12 36
13 36
14 36
23 36
24 36
25 36
28 36
29 36
30 36
35 36
36 36
37 36
40 36
41 36
42 36
46 36
47 36
48 36
55 36
56 36
57 36
58 36
59 36
70 36
71 36
72 36
78 36
79 36
3 37
4 37
5 37
12 37
13 37
14 37
23 37
24 37
28 37
29 37
35 37
36 37
37 37
40 37
41 37
42 37
45 37
46 37
47 37
57 37
58 37
59 37
70 37
71 37
72 37
78 37
79 37
3 38
4 38
5 38
12 38
13 38
14 38
22 38
23 38
24 38
28 38
29 38
35 38
36 38
37 38
40 38
41 38
42 38
45 38
46 38
47 38
57 38
58 38
59 38
69 38
70 38
71 38
78 38
79 38
3 39
4 39
5 39
11 39
12 39
13 39
14 39
22 39
23 39
24 39
27 39
28 39
29 39
35 39
36 39
37 39
40 39
41 39
44 39
45 39
46 39
57 39
58 39
59 39
69 39
70 39
71 39
78 39
79 39
80 39
3 40
4 40
5 40
11 40
12 40
13 40
22 40
23 40
24 40
27 40
28 40
29 40
35 40
36 40
37 40
40 40
41 40
42 40
43 40
44 40
45 40
57 40
58 40
59 40
69 40
70 40
71 40
78 40
79 40
80 40
3 41
4 41
11 41
12 41
13 41
21 41
22 41
23 41
27 41
28 41
29 41
35 41
36 41
37 41
40 41
41 41
42 41
43 41
56 41
57 41
58 41
59 41
60 41
69 41
70 41
71 41
78 41
79 41
80 41
3 42
4 42
10 42
11 42
12 42
21 42
22 42
23 42
28 42
29 42
35 42
36 42
40 42
41 42
42 42
56 42
57 42
58 42
59 42
60 42
69 42
70 42
71 42
78 42
79 42
80 42
3 43
4 43
10 43
11 43
12 43
21 43
22 43
28 43
29 43
35 43
36 43
40 43
41 43
42 43
55 43
56 43
57 43
59 43
60 43
61 43
69 43
70 43
71 43
78 43
79 43
80 43
3 44
4 44
9 44
10 44
11 44
20 44
21 44
22 44
28 44
29 44
30 44
34 44
35 44
36 44
40 44
41 44
42 44
55 44
56 44
57 44
60 44
61 44
69 44
70 44
77 44
78 44
79 44
3 45
4 45
9 45
10 45
11 45
20 45
21 45
29 45
30 45
34 45
35 45
40 45
41 45
42 45
55 45
56 45
57 45
60 45
61 45
69 45
70 45
77 45
78 45
79 45
2 46
3 46
4 46
9 46
10 46
11 46
19 46
20 46
21 46
29 46
30 46
34 46
35 46
40 46
41 46
42 46
54 46
55 46
56 46
60 46
61 46
62 46
68 46
69 46
70 46
77 46
78 46
79 46
2 47
3 47
9 47
10 47
11 47
18 47
19 47
20 47
21 47
30 47
31 47
32 47
33 47
34 47
35 47
40 47
41 47
42 47
54 47
55 47
56 47
61 47
62 47
68 47
69 47
70 47
76 47
77 47
78 47
2 48
3 48
9 48
10 48
11 48
17 48
18 48
19 48
20 48
21 48
22 48
30 48
31 48
32 48
33 48
34 48
41 48
42 48
53 48
54 48
55 48
61 48
62 48
68 48
69 48
70 48
76 48
77 48
78 48
2 49
3 49
4 49
9 49
10 49
11 49
16 49
17 49
18 49
21 49
22 49
31 49
32 49
33 49
41 49
42 49
52 49
53 49
54 49
61 49
62 49
63 49
67 49
68 49
69 49
75 49
76 49
77 49
2 50
3 50
4 50
9 50
10 50
11 50
14 50
15 50
16 50
17 50
21 50
22 50
23 50
32 50
33 50
34 50
41 50
42 50
51 50
52 50
53 50
54 50
62 50
63 50
67 50
68 50
69 50
75 50
76 50
77 50
3 51
4 51
10 51
11 51
12 51
13 51
14 51
15 51
22 51
23 51
24 51
32 51
33 51
34 51
41 51
42 51
51 51
52 51
53 51
62 51
63 51
67 51
68 51
74 51
75 51
76 51
3 52
4 52
5 52
10 52
11 52
12 52
13 52
14 52
23 52
24 52
33 52
34 52
41 52
42 52
43 52
51 52
52 52
62 52
63 52
64 52
65 52
66 52
67 52
68 52
74 52
75 52
76 52
4 53
5 53
6 53
11 53
12 53
13 53
23 53
24 53
33 53
34 53
41 53
42 53
43 53
51 53
52 53
63 53
64 53
65 53
66 53
73 53
74 53
75 53
4 54
5 54
6 54
11 54
12 54
13 54
23 54
24 54
25 54
33 54
34 54
42 54
43 54
51 54
52 54
63 54
64 54
65 54
66 54
73 54
74 54
75 54
4 55
5 55
6 55
12 55
13 55
14 55
23 55
24 55
25 55
33 55
34 55
42 55
43 55
51 55
52 55
63 55
64 55
65 55
73 55
74 55
75 55
4 56
5 56
6 56
12 56
13 56
14 56
23 56
24 56
25 56
33 56
34 56
35 56
42 56
43 56
51 56
52 56
53 56
63 56
64 56
65 56
73 56
74 56
75 56
4 57
5 57
6 57
12 57
13 57
14 57
15 57
23 57
24 57
25 57
33 57
34 57
35 57
42 57
43 57
44 57
51 57
52 57
53 57
63 57
64 57
65 57
73 57
74 57
75 57
4 58
5 58
6 58
13 58
14 58
15 58
23 58
24 58
25 58
33 58
34 58
35 58
42 58
43 58
44 58
52 58
53 58
54 58
63 58
64 58
65 58
73 58
74 58
75 58
4 59
5 59
14 59
15 59
16 59
22 59
23 59
24 59
33 59
34 59
43 59
44 59
53 59
54 59
63 59
64 59
65 59
74 59
75 59
4 60
5 60
15 60
16 60
22 60
23 60
33 60
34 60
43 60
44 60
53 60
54 60
55 60
64 60
65 60
66 60
74 60
75 60
76 60
4 61
5 61
15 61
16 61
17 61
21 61
22 61
23 61
33 61
34 61
44 61
45 61
54 61
55 61
65 61
66 61
67 61
75 61
76 61
77 61
4 62
5 62
15 62
16 62
17 62
21 62
22 62
23 62
32 62
33 62
34 62
44 62
45 62
54 62
55 62
56 62
66 62
67 62
68 62
75 62
76 62
77 62
4 63
5 63
16 63
17 63
21 63
22 63
32 63
33 63
34 63
44 63
45 63
54 63
55 63
56 63
66 63
67 63
68 63
76 63
77 63
78 63
4 64
5 64
16 64
17 64
18 64
20 64
21 64
22 64
32 64
33 64
34 64
44 64
45 64
46 64
55 64
56 64
67 64
68 64
69 64
77 64
78 64
4 65
5 65
16 65
17 65
18 65
20 65
21 65
22 65
31 65
32 65
33 65
44 65
45 65
46 65
55 65
56 65
57 65
67 65
68 65
69 65
77 65
78 65
79 65
4 66
5 66
6 66
16 66
17 66
18 66
20 66
21 66
30 66
31 66
32 66
44 66
45 66
46 66
56 66
57 66
67 66
68 66
69 66
78 66
79 66
80 66
4 67
5 67
6 67
16 67
17 67
18 67
19 67
20 67
21 67
30 67
31 67
32 67
45 67
46 67
56 67
57 67
58 67
67 67
68 67
69 67
79 67
80 67
4 68
5 68
6 68
16 68
17 68
18 68
19 68
20 68
29 68
30 68
31 68
45 68
46 68
57 68
58 68
59 68
67 68
68 68
69 68
79 68
80 68
81 68
5 69
6 69
16 69
17 69
18 69
19 69
28 69
29 69
30 69
31 69
44 69
45 69
46 69
47 69
57 69
58 69
59 69
60 69
67 69
68 69
69 69
79 69
80 69
81 69
5 70
6 70
15 70
16 70
17 70
18 70
28 70
29 70
30 70
44 70
45 70
46 70
58 70
59 70
60 70
67 70
68 70
79 70
80 70
81 70
5 71
6 71
15 71
16 71
17 71
27 71
28 71
29 71
30 71
44 71
45 71
46 71
59 71
60 71
61 71
67 71
68 71
79 71
80 71
5 72
6 72
14 72
15 72
16 72
26 72
27 72
29 72
30 72
44 72
45 72
46 72
59 72
60 72
61 72
66 72
67 72
68 72
78 72
79 72
80 72
4 73
5 73
6 73
13 73
14 73
15 73
25 73
26 73
27 73
29 73
30 73
31 73
43 73
44 73
45 73
59 73
60 73
61 73
66 73
67 73
68 73
78 73
79 73
4 74
5 74
6 74
12 74
13 74
14 74
15 74
25 74
26 74
27 74
30 74
31 74
43 74
44 74
45 74
59 74
60 74
61 74
66 74
67 74
68 74
77 74
78 74
79 74
4 75
5 75
11 75
12 75
13 75
14 75
24 75
25 75
26 75
30 75
31 75
32 75
43 75
44 75
45 75
59 75
60 75
61 75
66 75
67 75
76 75
77 75
78 75
3 76
4 76
5 76
11 76
12 76
13 76
23 76
24 76
25 76
31 76
32 76
43 76
44 76
45 76
59 76
60 76
61 76
66 76
67 76
76 76
77 76
78 76
3 77
4 77
11 77
12 77
13 77
22 77
23 77
24 77
31 77
32 77
33 77
43 77
44 77
45 77
59 77
60 77
61 77
65 77
66 77
67 77
75 77
76 77
77 77
2 78
3 78
4 78
11 78
12 78
13 78
22 78
23 78
32 78
33 78
42 78
43 78
44 78
45 78
58 78
59 78
60 78
65 78
66 78
67 78
75 78
76 78
77 78
1 79
2 79
3 79
11 79
12 79
13 79
21 79
22 79
23 79
32 79
33 79
34 79
42 79
43 79
44 79
57 79
58 79
59 79
60 79
65 79
66 79
75 79
76 79
77 79
1 80
2 80
3 80
11 80
12 80
13 80
14 80
20 80
21 80
22 80
32 80
33 80
34 80
43 80
44 80
57 80
58 80
59 80
65 80
66 80
75 80
76 80
77 80
1 81
2 81
12 81
13 81
14 81
15 81
19 81
20 81
32 81
33 81
34 81
43 81
44 81
45 81
56 81
57 81
58 81
64 81
65 81
66 81
75 81
76 81
77 81
1 82
2 82
13 82
14 82
15 82
16 82
17 82
18 82
19 82
32 82
33 82
34 82
44 82
45 82
46 82
55 82
56 82
57 82
58 82
64 82
65 82
75 82
76 82
77 82
14 83
15 83
16 83
17 83
18 83
32 83
33 83
34 83
45 83
46 83
47 83
55 83
56 83
57 83
64 83
65 83
75 83
76 83
77 83
15 84
16 84
17 84
32 84
33 84
34 84
46 84
47 84
48 84
54 84
55 84
56 84
57 84
64 84
65 84
75 84
76 84
    ];

A = sparse(IJ(:,1),IJ(:,2),1);
A = full(A);