% Program to carry out facetted imaging of a given visibility set.
% The gridding and FFT is from fft_imager_sjw.m

% Imaging of sampled visibility using FFT, after gridding the sampled 
% visibility onto a uniform grid.
% SJW, with additions by Peeyush

% Parameters:
%  acc = observed visibility at coords. specified by u,v below, full matrix for all unflagged antennas
%  u   = u-coord of measured visibility in local coordinates, in m
%  v   = v-coord of measured visibility in local coordinates, in m
%  duv = uv-coord. resolution of desired uv grid, in m
%  Nuv = No. of points in grid, dimensionless
% uvsize = Final size of grid, dimensionless, so that actual final size in m =
% uvsize * duv
% freq = Freq. of observation, in Hz
% nfacets = Number of facets to divide the desired FoV into

% pep/26Mar12

function [skymap, vispad, l, m] = fft_imager_facet (acc, u, v, duv, Nuv, uvsize, freq, nfacets)
   %% Currently supporting only4,9,16, 25 ...  facets
   if (sqrt (nfacets) - floor (sqrt(nfacets)) > 0)
       disp (['Cant create non-square number of facets!!']);
       return;
   end
   
   if uvsize < Nuv
       disp (['Cant have padded image smaller than FFT size!']);
       uvsize = Nuv;
   end
   % Determine the FoV of a single image, based on specified parameters,
   % in l,m units
   dlmin = (299792458/(uvsize*duv*freq)); % best l,m resolution
   lmax  = dlmin*uvsize/2;                % Max extent of l,m
   
   %% Determine the facet centers in local l,m coordinates, assuming l=0,m=0
   % is at the phase center for a zenith pointing (ie, at center of lm
   % grid), and FoV based on specified parameters
   nrows = sqrt (nfacets); ncols = nrows; % NOTE: square facet grids only!
   facet_cen = zeros (2, nfacets); % l,m coord. of facet centers
   for rind = 1:nrows
       l_row = -lmax + 1/nrows*(2*(rind-1) + 1); % l coord. of all facets in this row
       for cind = 1:ncols
           facet_cen (1, (rind-1)*nrows + cind) = l_row;
           facet_cen (2, (rind-1)*nrows + cind) = -lmax + 1/ncols*(2*(cind-1) + 1);   % l coord. of facet       
       end
   end
   
   % Determine l,m coords of facets, such that image resolution of a facet
   % is identical to that of a single, all sky image
   lfacet = lmax/nrows;        % this is half the available l-range, per facet
   pix2facet = 2*lfacet/dlmin; % Total number of pixels per facet
   duv_facet = (299792458/(freq*pix2facet*dlmin));
   l = zeros (pix2facet, nfacets);
   m = l;
   for ind = 1:nfacets
       lcen = facet_cen (1, ind); mcen = facet_cen (2, ind);
       l(:, ind) = linspace (lcen-lfacet, lcen+lfacet, pix2facet);
       m(:, ind) = linspace (mcen-lfacet, mcen+lfacet, pix2facet);
   end
   
   disp (['(l,m) Coordinates of facet centers: ']);
   for ind = 1:nrows*ncols
       disp (['(' num2str(facet_cen (1, ind), '%4.2f') ',' num2str(facet_cen(2,ind), '%4.2f') ')']);
   end
    
   % skymap = zeros (Nuv, Nuv, nfacets); 
   % vis = zeros (Nuv, Nuv, nfacets);
   
   skymap = zeros (pix2facet, pix2facet, nfacets);
   vis = zeros (pix2facet, pix2facet, nfacets);
   
   %% Gridding UV data uniformly, for every facet
   for ind = 1:2
    for idx = 1:length(u(:)) 
        ampl = abs(acc(idx));
        phasor = acc(idx) / ampl;
        uidx = u(idx) / duv_facet + pix2facet / 2;
        uidxl = floor(uidx);
        uidxh = ceil(uidx);
        dul = abs(uidx - uidxl);
        duh = abs(uidx - uidxh);
        sul = duh * ampl;
        suh = dul * ampl;

        vidx = v(idx) / duv_facet + pix2facet / 2;
        vidxl = floor(vidx);
        vidxh = ceil(vidx);
        dvl = abs(vidx - vidxl);
        dvh = abs(vidx - vidxh);
        sull = dvh * sul;
        sulh = dvl * sul;
        suhl = dvh * suh;
        suhh = dvl * suh;

        vis(uidxl, vidxl, ind) = vis(uidxl, vidxl, ind) + sull * phasor;
        vis(uidxl, vidxh, ind) = vis(uidxl, vidxh, ind) + sulh * phasor;
        vis(uidxh, vidxl, ind) = vis(uidxh, vidxl, ind) + suhl * phasor;
        vis(uidxh, vidxh, ind) = vis(uidxh, vidxh, ind) + suhh * phasor;

        %W(uidx, vidx) = W(uidx, vidx) + 1;
    end
    vis (~isfinite(vis(:,:,ind))) = 0;
    
    % Image the gridded facets
    skymap (:,:,ind)  = fftshift (fft2 (vis(:,:,ind)));
   end

%% zero padding to desired (u,v)-size
% N = size(vis, 1);
% N1 = floor((uvsize - N) / 2);
% N2 = ceil((uvsize + 1 - N) / 2) - 1;
% vispad = [zeros(N1, uvsize); ...
%           zeros(N, N1), vis, zeros(N, N2); ...
%           zeros(N2, uvsize)];
% vispad(~isfinite(vispad)) = 0;
% 
% % compute image
% ac = zeros(size(vispad));
% ac(256, 256) = 100;
% vispad = vispad + ac;
% vispad = conj(flipud(fliplr(fftshift(vispad))));

% %% Generate images per facet
% for ind = 1:nrows*ncols
%     % Rephase array in the direction of a facet
%     ph = exp (-i*2*pi*(facet_cen(1, ind)*pix2lm*uloc/Nuv + facet_cen(2, ind)*pix2lm*vloc/Nuv);
%     vispad_shift = vispad .* ph;
%     skymap (:,:,ind)  = fftshift (fft2 (vispad_shift));
% end