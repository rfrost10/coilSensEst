function [csm] = ismrm_estimate_csm_walsh_3d(img, smoothing)
%
%   [csm] = ismrm_estimate_csm_walsh(img)
%
%   Estimates relative coil sensitivity maps from a set of coil images
%   using the eigenvector method described by Walsh et al. (Magn Reson Med
%   2000;43:682-90.)
%
%   INPUT:
%     - img     [x,y,z,coil]   : Coil images
%     - smooth  scalar       : Smoothing block size (defaults to 5)
%
%   OUTPUT:
%
%     - csm     [x,y,z,coil    : Relative coil sensitivity maps
%
%
%   Code is based on an original implementation by Peter Kellman, NHLBI,
%   NIH
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%
%   June 2020: RobF: modified to work for 3D images

if nargin < 2,
    smoothing = 5;
end

ncoils = size(img,4);
nslices = size(img,3);

% normalize by root sum of squares magnitude
mag = sqrt(sum(img .* conj(img),4));
s_raw=img ./ repmat(mag + eps,[1 1 1 ncoils]); clear mag; 


% compute sample correlation estimates at each pixel location
Rs=ismrm_correlation_matrix(s_raw);


% apply spatial smoothing to sample correlation estimates (NxN convolution)
if smoothing>1,
	h_smooth=ones(smoothing)/(smoothing^2); % uniform smoothing kernel
    for s=1:nslices
        disp(s)
	for m=1:ncoils
		for n=1:ncoils
		    Rs(:,:,s,m,n)=conv2(Rs(:,:,s,m,n),h_smooth,'same');
		end
    end
    end
end


% compute dominant eigenvectors of sample correlation matrices
[csm,lambda]=ismrm_eig_power(Rs); % using power method


return 


%Utility functions provided by Peter Kellman, NIH.
function [Rs]=ismrm_correlation_matrix(s)
% function [Rs]=correlation_matrix(s);
%
% function correlation_matrix calculates the sample correlation matrix (Rs) for
% each pixel of a multi-coil image s(y,x,coil)
%
% input:
%    s   complex multi-coil image s(y,x,coil)
% output:
%    Rs  complex sample correlation matrices, Rs(y,x,coil,coil)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

[rows,cols,slices,ncoils]=size(s);
Rs=zeros(rows,cols,slices,ncoils,ncoils); % initialize sample correlation matrix to zero
for i=1:ncoils
    for j=1:i-1
		Rs(:,:,:,i,j)=s(:,:,:,i).*conj(s(:,:,:,j));
		Rs(:,:,:,j,i)=conj(Rs(:,:,:,i,j)); % using conjugate symmetry of Rs
    end
	Rs(:,:,:,i,i)=s(:,:,:,i).*conj(s(:,:,:,i));
end

return

function [v,d]=ismrm_eig_power(R);
% function [v,d]=eig_power(R);
%
% vectorized method for calculating the dominant eigenvector based on
% power method. Input, R, is an image of sample correlation matrices
% where: R(y,x,z,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%
% v is the dominant eigenvector
% d is the dominant (maximum) eigenvalue

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

rows=size(R,1);cols=size(R,2);slices=size(R,3);ncoils=size(R,4);
N_iterations=2;
v=ones(rows,cols,slices,ncoils); % initialize e.v.

d=zeros(rows,cols,slices);
for i=1:N_iterations
    v=squeeze(sum(R.*repmat(v,[1 1 1 1 ncoils]),4));
	d=ismrm_rss(v);
    d( d <= eps) = eps;
	v=v./repmat(d,[1 1 1 ncoils]);
end

p1=angle(conj(v(:,:,:,1)));
% (optionally) normalize output to coil 1 phase
v=v.*repmat(exp(sqrt(-1)*p1),[1 1 1 ncoils]);
v=conj(v);

return