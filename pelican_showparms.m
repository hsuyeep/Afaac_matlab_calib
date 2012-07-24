% Script to generate images projected in RA/Dec plane, and to show
% gainsolutions and other parameters.
% Operates on the output of pelican_pipesim_storevars.m
% pep/19Apr12

%% 
dirname = './data/';
% [~, lsout] = dos(['ls -1 ' dirname '*clean.txt']);
[~, lsout] = dos(['ls -1 ' dirname '48233*.mat']);

datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
disp (['Found ' num2str(nfiles) ' files']);
Nelem = 288; % NOTE Hardcoded

% Determine the l,m coordinates, based on parameters chosen
% NOTE: imager used was fft_imager.m, with parameters below:
duv = 2; Nuv = 1000; uvsize = 1024; freq = 59756469; c = 299792458;
lambda = c / freq;
ltick_res = lambda / (uvsize * duv);
ltick_max = uvsize * ltick_res / 2;
ltick = -ltick_max:ltick_res:ltick_max-1e-3;
mask = NaN(length(ltick));
mask(meshgrid(ltick).^2 + meshgrid(ltick).'.^2 < 1) = 1;

figure
imagesc(ltick, ltick, real(calmap) .* mask);
set(gca, 'FontSize', 16);
title('calibrated sky map, extended emission removed');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow xl \rightarrow West');
set(colorbar, 'FontSize', 16);
axis([-1 1 -1 1]);

return ;
lfft = linspace (-1, 1, uvpad);
mfft = lfft;

sky_tobs = zeros (nfiles);
sky_l = zeros (z(2), nfiles);
sky_m = zeros (z(2), nfiles);


fnum = 1;
%% 
for fnum = 1:nfiles
    disp (['Processing file ' num2str(fnum) ' of ' num2str(nfiles)]);
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);
    
    %% convert images to RA/dec coordinate system
    % First convert l.m coordinates to RA/dec at the time of observation
    [alpha, delta] = lmtoradec (lfft, mfft, JulianDay (tobs));
    sel = ~isnan (alpha(:));
    
    % Then interpolate the image to create a model surface
    radecimage = TriScatteredInterp (alpha (sel), delta (sel), abs(tmp_image (sel)));
    radecimage_subAteam = TriScatteredInterp (alpha(sel), delta(sel), abs (tmp_image(sel)));
    % Create the regularly sampled RA/dec plane
    [ragrid, decgrid] = meshgrid (linspace (0,2*pi, uvpad), linspace (-pi/2,pi/2, uvpad));
    
    % generate samples from model skyimage
    sky_radecimage(:,:,fnum) = radecimage (ragrid, decgrid);
    
    % Do the same for A-team subtracted images
    [tmp_image ~] = fft_imager_sjw (accsubAteam, uloc (:), vloc (:), duv ,Nuv, uvpad);  
    sky_subAteam (:, :, fnum) = tmp_image;
    
    % Then interpolate the image to create a model surface
    radecimage_subAteam = TriScatteredInterp (alpha(sel), delta(sel), abs (tmp_image(sel)));
    
    % generate samples from model skyimage
    sky_radecimage_subAteam(:,:,fnum) = radecimage_subAteam (ragrid, decgrid);
    
    sky_alpha (:,:,fnum) = alpha;
    sky_delta (:,:,fnum) = delta;
    % sky_subAteam (:, :, fnum) = acm2skyimage(accsubAteam, poslocal(:, 1), poslocal(:, 2), freq, l, m);
  
    
    
%     imagesc(l, m, skymapcal1 .* mask);
%     set(gca, 'FontSize', 16);
%     % title('extended emission, Sun and A-team removed');
%     axis equal
%     axis tight
%     xlabel('South \leftarrow m \rightarrow North');
%     ylabel('East \leftarrow l \rightarrow West');
%     set(colorbar, 'FontSize', 16);
end

%%
outfilename = 'sky_radec_images.mat';
save (outfilename, 'sky_radecimage', 'sky_tobs', 'ragrid', 'sky_m');

% outfilename = 'sky_subAteam_images.mat';
%save (outfilename, 'sky_subAteam', 'sky_tobs', 'sky_l', 'sky_m');
