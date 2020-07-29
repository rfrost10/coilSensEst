function [recon,cmap,wfull]=adapt_array_3d(yn,rn,norm);
% function [recon,cmap,wfull]=adapt_array_3d(yn,rn,norm)
% Reconstruction of array data and computation of coil sensitivities based 
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%-------------------------------------------------------------------------
%	Input:
%	yn: array data to be combined [ny, nx, nz, nc]. 
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	recon: reconstructed image [ny, nx, nz].
%	cmap: estimated coil sensitivity maps [ny, nx, nz, nc].
%--------------------------------------------------------------------------
% Ricardo Otazo
% CBI, New York University
%--------------------------------------------------------------------------
%
% 20200513: adapted by Rob Frost for 3D data, starting from adapt_array_2d

yn=permute(yn,[4,1,2,3]);
[nc,ny,nx,nz]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

% find coil with maximum intensity for phase correction
[mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[3 2 4 1])))));   

bs1=8;  %x-block size
bs2=8;  %y-block size
bs3=8;  %z-block size
st=4;   %increase to set interpolation step size

% RobF why did x not have round() here?
wsmall=zeros(nc,round(ny./st), round(nx./st), round(nz./st)); 
cmapsmall=zeros(nc,round(ny./st), round(nx./st), round(nz./st));

tic
for z=st:st:nz
for x=st:st:nx
for y=st:st:ny
    %Collect block for calculation of blockwise values
    ymin1=max([y-bs1./2 1]);                   
    xmin1=max([x-bs2./2 1]);
    zmin1=max([z-bs3./2 1]);
    % Cropping edges
    ymax1=min([y+bs1./2 ny]);                 
    xmax1=min([x+bs2./2 nx]);                  
    zmax1=min([z+bs3./2 nz]);

    ly1=length(ymin1:ymax1);
    lx1=length(xmin1:xmax1);
    lz1=length(zmin1:zmax1);
    m1=reshape(yn(:,ymin1:ymax1,xmin1:xmax1,zmin1:zmax1),nc,lx1*ly1*lz1);
      
    m=m1*m1'; %signal covariance
      
    % eignevector with max eigenvalue for optimal combination
    [e,v]=eig(inv(rn)*m);                    
                                               
    v=diag(v);
    [mv,ind]=max(v);
      
    mf=e(:,ind);                      
    mf=mf/(mf'*inv(rn)*mf);               
    normmf=e(:,ind);
    
    % Phase correction based on coil with max intensity
    mf=mf.*exp(-j*angle(mf(maxcoil)));        
    normmf=normmf.*exp(-j*angle(normmf(maxcoil)));

    wsmall(:,y./st,x./st,z./st)=mf;
    cmapsmall(:,y./st,x./st,z./st)=normmf;
end
end
end
disp(['creating weights took ' num2str(toc) ' sec'])

recon=zeros(ny,nx,nz);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude 
% pixels between +1 and -1 pixels.
tic
for i=1:nc
    % RobF: can't use 'bilinear' as in adapt_array_2d.m --> use 'cubic'
wfull(i,:,:,:)=conj(imresize3(squeeze(abs(wsmall(i,:,:,:))),[ny nx nz],'cubic') ...        
.*exp(j.*imresize3(angle(squeeze(wsmall(i,:,:,:))),[ny nx nz],'nearest')));
cmap(i,:,:,:)=imresize3(squeeze(abs(cmapsmall(i,:,:,:))),[ny nx nz],'cubic') ...
.* exp(j.*imresize3(squeeze(angle(cmapsmall(i,:,:,:))),[ny nx nz],'nearest'));
end
disp(['interpolating weights to full res. took ' num2str(toc) ' sec'])

recon=squeeze(sum(wfull.*yn));   %Combine coil signals. 
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);
wfull = permute(wfull,[2,3,4,1]);