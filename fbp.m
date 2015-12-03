function x = fbp(A,b,theta,filter)
%FBP  Filtered Back Projection that uses the sparse system matrix
%
% x = fbp(A,b,theta)
% x = fbp(A,b,theta,filter)
%
% Input: A = system matrix,
%        b = right-hand side (vectorized sinogram),
%        theta = vector of projection angles,
%        filter = see iradon for details.
%
% This function in quite similar to Matlab's iradon, but it is designed to
% conform with AIR Tools. The filter is identical to that in iradon.  The
% backprojection in performed by multiplication with the transpose of A,
% and the sinogram is represented as the vector b.

% Per Christian Hansen, DTU Compute, December 12, 2014.

% Input check.
if nargin < 3, error('Input A, b, theta must me specified'); end

% Default filter
if nargin < 4, filter = 'Ram-Lak'; end
filter = lower(filter);  % Turn all characters to lower case.

if strcmp(filter,'none'), x = A'*b; return, end

% Reshape right-hand side b to a sinogram.
m = size(A,1);       % Number of rows in A.
na = length(theta);  % Number of angles.
nr = m/na;           % Number of rays.
S = reshape(b,nr,na);

% Filter the sinogram.
d = 1;
S = filterProjections(S, filter, d);

% Bckprojection.
x = A'*S(:)*pi/(2*length(theta));

%========== These subfunctions are from iradon ========================

function [p,H] = filterProjections(p_in, filter, d)

p = p_in;

% Design the filter
len = size(p,1);
H = designFilter(filter, len, d);

if strcmpi(filter, 'none')
    return;
end

p(length(H),1)=0;  % Zero pad projections

% In the code below, I continuously reuse the array p so as to
% save memory.  This makes it harder to read, but the comments
% explain what is going on.

p = fft(p);    % p holds fft of projections

for i = 1:size(p,2)
    p(:,i) = p(:,i).*H; % frequency domain filtering
end

p = real(ifft(p));     % p is the filtered projections
p(len+1:end,:) = [];   % Truncate the filtered projections

function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections

order = max(64,2^nextpow2(2*len));

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter