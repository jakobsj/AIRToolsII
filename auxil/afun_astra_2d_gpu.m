function y = afun_astra_2d_gpu(x,transp_flag,proj_geom,vol_geom)
%AFUN_ASTRA_2D_GPU  Wrap ASTRA projectors into a function handle
%
%   y = afun_astra_2d_gpu(x,transp_flag,proj_geom,vol_geom)
%
% This function takes an ASTRA projection geometry proj_geom and volume
% geometry vol_geom specifying the configuration of tomographic forward and 
% backprojection, and wraps them into a function in the format expected by
% AIR Tools functions. afun_astra_2d_gpu uses the GPU-mode of ASTRA and is
% compatible with all 2D cases implemented in ASTRA. Please see the ASTRA
% documentation for further details, http://www.astra-toolbox.com/.
% 
% Depending on the second input, either the system matrix or its transpose
% is multiplied onto the first input vector, or the size of the system
% matrix is returned.
%
% Given ASTRA proj_geom and vol_geom, this function wraps them in an
% anonymous function
%    myfun = @(XX,TT) afun_astra_2d_gpu(XX,TT,proj_geom,vol_geom);
% after which myfun implicitly can apply A or A' or return the size of A,
% implemented using ASTRA on the GPU, by calling:
%    y = myfun(x,'notransp');
%    z = myfun(y,'transp');
%    s = myfun([],'size');
%
% Input:
%   x            Vector on which to apply matrix multiplication from the
%                left by either A or A'; x must be a column vector with
%                length matching the relevant dimension of A.
%   transp_flag  String to indicate whether to apply multiplication by A
%                ('notransp') or A' ('transp'), or return the size of A
%                ('size'). If set to 'size', the input x is ignored.
%   proj_geom    ASTRA projection geometry.
%   vol_geom     ASTRA volume geometry.
%
% Output:
%   y            If transp_flag is 'notransp' then y is A*x; if 'transp'
%                then y is A'*x; if 'size' then y is size(A).
%
% See also: afun_matrix, demo_astra_2d.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.
% Tested using ASTRA 1.8 on Linux.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

switch transp_flag
    case 'size'
        % Extract dimension of operator from proj_geom and vol_geom
        y = [proj_geom.DetectorCount*length(proj_geom.ProjectionAngles),...
             vol_geom.GridRowCount*vol_geom.GridColCount];
    case 'notransp'
        % x is vectorized image to be reshaped. Different geometry
        % conventions used in AIR Tools and ASTRA, so need to rotate by 180
        % degrees, and transpose and flip sinogram created by ASTRA to
        % match with AIR Tools. Vectorize before returning y as output.
        X = rot90(reshape(x,vol_geom.GridRowCount,vol_geom.GridColCount),2);
        [sid,y] = astra_create_sino_cuda(X,proj_geom,vol_geom);
        astra_mex_data2d('delete', sid);
        y = flipud(y');
        y = y(:);
    case 'transp'
        % x is vectorized sinogram to be reshaped. As above, need to flip
        % and transpose sinogram before backprojecting, and then rotate
        % image by 180 degrees to match with AIR Tools conventions. 
        % Vectorize before returning y as output.
        x = reshape(x,proj_geom.DetectorCount,...
                      length(proj_geom.ProjectionAngles));
        x = flipud(x)';
        y = astra_create_backprojection_cuda(x,proj_geom,vol_geom);
        y = rot90(y,2);
        y = y(:);
    otherwise
        error('transp_flag must be ''size'', ''notransp'' or ''transp''')
end