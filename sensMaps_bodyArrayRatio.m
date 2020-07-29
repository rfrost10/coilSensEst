function [maps,maps_masked] = sensMaps_bodyArrayRatio(im_mc,im_bc,mask)
% [poly_img] = sensMaps_bodyArrayRatio(im_mc,im_bc,mask)
%
% fit sensitivity maps to ratio of multi- and single-coil images
%
% INPUTS:
%       im_mc - multi-channel images [nx,ny,nz,nc]
%       im_bc - body coil (single-channel) images [nx,ny,nz]
%
% OUTPUTS:
%       maps - sensitivity maps
%       maps_masked - masked sensitivity maps

if length(size(im_mc)) ~= 4
    size(im_mc)
    error('multi-channel data should be 4D array!')
end

% dimensions
[nx,ny,nz,nc] = size(im_mc);

if any(size(im_bc) ~= [nx,ny,nz])
    disp(size(im_mc)); disp(size(im_bc))
    error('single-channel data should match the size of multi-channel images')
end
if any(size(mask) ~= [nx,ny,nz])
    disp(size(im_mc)); disp(size(mask))
    error('mask should match the size of multi-channel images')
end

[xx, yy, zz] = ndgrid(1:nx, 1:ny, 1:nz);
xx(~mask) = [];
yy(~mask) = [];
zz(~mask) = [];

% scale up the body coil image
im_mc_sos = sqrt(sum(abs(im_mc).^2, 4));
scaleFac = mean(im_mc_sos(mask>0)) / mean(abs(im_bc(mask>0)));
im_bc = im_bc * scaleFac;

mapsRatio = im_mc ./ repmat(im_bc, [1 1 1 nc]);

[xx_eval, yy_eval, zz_eval] = ndgrid(1:nx, 1:ny,1:nz);
% matrices for magnitude (real) fit
depvar_re = zeros(nx, ny, nz);
poly_re = zeros(nx*ny*nz,nc);

% matrices for phase (imaginary) fit
depvar_im=zeros(nx, ny, nz);
poly_im = zeros(nx*ny*nz,nc);

% Choose polynomial order
order=7;

tic
for cha=1:nc
    % Polynomial fit for real parts
    depvar_re = real(mapsRatio(:,:,:,cha));
    depvar_re(~mask) = [];
    indepvar_real = polyfitn([xx(:), yy(:), zz(:)], depvar_re(:), order);
    tmp = polyvaln(indepvar_real, [xx_eval(:), yy_eval(:),zz_eval(:)]);
    poly_re(:,cha) = tmp;
    
    % Polynomial fit for imaginary parts
    depvar_im = imag(mapsRatio(:,:,:,cha));
    depvar_im(~mask) = [];
    indepvar_imag = polyfitn([xx(:), yy(:), zz(:)], depvar_im(:), order);
    tmp = polyvaln(indepvar_imag, [xx_eval(:), yy_eval(:),zz_eval(:)]);
    poly_im(:,cha) = tmp;        
end
mapTime = toc;
disp(['sensitivity fitting time: ' num2str(mapTime) ' sec' ])

% rewrite to complex numbers from magnitude and phase
poly_re = reshape(poly_re, [nx ny nz nc]);
poly_im = reshape(poly_im, [nx ny nz nc]);
maps = real(poly_re) + 1i*real(poly_im);

maps_masked = maps .* repmat(abs(mask),[1 1 1 nc]);