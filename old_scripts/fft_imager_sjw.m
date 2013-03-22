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

function [skymap, vispad, l, m] = fft_imager_sjw(acc, u, v, duv, Nuv, uvsize, freq)

% Determine size of visibility grid, based on duv
Nuv = ceil (max (max (u), max (v))/duv*2) + 2; % +2 for accomodating both ends of u/v range
disp (['Nuv = ' num2str(Nuv) 'max (u,v)= ' num2str(max(max(u))) num2str(max(max(v)))]);

% create object for interpolation
vis = zeros(Nuv);

%W = zeros(Nuv);
for idx = 1:length(u(:)) 
    ampl = abs(acc(idx));
    phasor = acc(idx) / ampl;
    uidx = u(idx) / duv + Nuv / 2;
    uidxl = floor(uidx);
    uidxh = ceil(uidx);
    dul = abs(uidx - uidxl);
    duh = abs(uidx - uidxh);
    sul = duh * ampl;
    suh = dul * ampl;
    
    vidx = v(idx) / duv + Nuv / 2;
    vidxl = floor(vidx);
    vidxh = ceil(vidx);
    dvl = abs(vidx - vidxl);
    dvh = abs(vidx - vidxh);
    sull = dvh * sul;
    sulh = dvl * sul;
    suhl = dvh * suh;
    suhh = dvl * suh;
    
    vis(uidxl, vidxl) = vis(uidxl, vidxl) + sull * phasor;
    vis(uidxl, vidxh) = vis(uidxl, vidxh) + sulh * phasor;
    vis(uidxh, vidxl) = vis(uidxh, vidxl) + suhl * phasor;
    vis(uidxh, vidxh) = vis(uidxh, vidxh) + suhh * phasor;
    
    %W(uidx, vidx) = W(uidx, vidx) + 1;
end

% zero padding to desired (u,v)-size
if uvsize > Nuv
    N = size(vis, 1);
    N1 = floor((uvsize - N) / 2);
    N2 = ceil((uvsize + 1 - N) / 2) - 1;
    vispad = [zeros(N1, uvsize); ...
              zeros(N, N1), vis, zeros(N, N2); ...
              zeros(N2, uvsize)];
else
    vispad = vis; % No padding requested
end
vispad(~isfinite(vispad)) = 0;

% Create l,m axis corresponding to choices of duv
dl = (299792458/(freq * uvsize * duv)); % dimensionless, in dir. cos. units
lmax = dl * uvsize / 2;
l = [-lmax:dl:lmax-1e-3]; 
m = l;  % Identical resolution and extent along m-axis

% compute image
ac = zeros(size(vispad));
ac(256, 256) = 100;
%vispad = vispad + ac;
vispad = conj(flipud(fliplr(fftshift(vispad))));
skymap = fftshift(fft2(vispad));