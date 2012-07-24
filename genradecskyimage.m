% Script to generate skyimages projected on the RA/DEC space from
% calibrated visibilities. The data used is the 10min sampling of the 3hr
% dataset from the first AFAAC offline correlations.
% pep/14Mar12

%% 
dirname = '../cookdat/per_snapshot_cal_sol/';
% [~, lsout] = dos(['ls -1 ' dirname '*clean.txt']);
[~, lsout] = dos(['ls -1 ' dirname '*clean_var.mat']);

datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
disp (['Found ' num2str(nfiles) ' files']);
Nelem = 288; % NOTE Hardcoded

% Needed do be done only once: load ITRF positions of antennas
fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
posITRF = zeros(Nelem, 3);
posITRF(:, 1) = posfiledata{3};
posITRF(:, 2) = posfiledata{4};
posITRF(:, 3) = posfiledata{5};

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF * rotmat;
uloc = meshgrid (poslocal (:,1)) - meshgrid (poslocal(:,1)).';
vloc = meshgrid (poslocal (:,2)) - meshgrid (poslocal(:,2)).';
% duv = 600/256;
% Nuv = 256;
% uvpad = 256; % Final FFT image will have this dimension
duv = 2;
Nuv = 1000;
uvpad = 1024;

% create 3D matrix to store images in time order
ccm_fname = datafiles {1}{1};
load (ccm_fname);
 
z = size (l);
% NOTE: if non FFT imaging is required, set size of sky_cal2 to zeros
% (z(2), z(2),nfiles);
sky_cal2 = zeros (uvpad, uvpad, nfiles);
sky_subAteam = zeros (uvpad, uvpad, nfiles);
sky_radecimage = zeros (uvpad, uvpad, nfiles);
sky_radecimage_subAteam = zeros (uvpad, uvpad, nfiles);

lfft = linspace (-1, 1, uvpad);
mfft = lfft;
alpha = zeros (uvpad, nfiles);
delta = zeros (uvpad, nfiles);


sky_tobs = zeros (nfiles);
sky_l = zeros (z(2), nfiles);
sky_m = zeros (z(2), nfiles);


fnum = 1;
%% 

outfilename = 'sky_radec_images_roysoc.mat';
try

for fnum = 1:nfiles
    disp (['Processing file ' num2str(fnum) ' of ' num2str(nfiles)]);
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);
    
    mask = NaN(length(l));
    mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;

    sky_tobs (fnum) = tobs; % NOTE: Time of observation is secs in Matlab datenum units! (secs from 01Jan0000) 
    sky_l (:, fnum)  = l;
    sky_m (:, fnum) = m;
    % Generate the skyimage with a DFT imager in l,m space
    % sky_cal2(:, :, fnum) = acm2skyimage(acccal2, poslocal(:, 1), poslocal(:, 2), freq, l, m);
    
% -------- calibrated image -----
    % Generate a skyimage with an FFT imager in l, m space
    [tmp_image ~] = fft_imager (acccal2(:), uloc (:), vloc (:), duv ,Nuv, uvpad);  
    sky_cal2 (:, :, fnum) = tmp_image;
    
    %% convert images to RA/dec coordinate system
    % First convert l.m coordinates to RA/dec at the time of observation
    [alpha, delta] = lmtoradec (lfft, mfft, JulianDay (tobs));
    sel = ~isnan (alpha(:));
    
    % Then interpolate the image to create a model surface
    radecimage = TriScatteredInterp (alpha (sel), delta (sel), abs(tmp_image (sel)));
    % Create the regularly sampled RA/dec plane
    [ragrid, decgrid] = meshgrid (linspace (0,2*pi, uvpad), linspace (-pi/2,pi/2, uvpad));
    
    % generate samples from model skyimage
    sky_radecimage(:,:,fnum) = radecimage (ragrid, decgrid);

% -------- calibrated image, A-team subtracted -----
    
    % Do the same for A-team subtracted images
    [tmp_image ~] = fft_imager (accsubAteam, uloc (:), vloc (:), duv ,Nuv, uvpad);  
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
catch err
        disp ('Error encountered! Saving variables to disk..');
        save (outfilename, 'sky_radecimage', 'sky_cal2', 'sky_alpha', 'sky_delta', 'sky_radecimage_subAteam', 'sky_subAteam', 'sky_tobs', 'ragrid', 'sky_m');
  
end 
%%
save (outfilename, 'sky_radecimage', 'sky_cal2', 'sky_alpha', 'sky_delta', 'sky_radecimage_subAteam', 'sky_subAteam', 'sky_tobs', 'ragrid', 'sky_m');

% outfilename = 'sky_subAteam_images.mat';
%save (outfilename, 'sky_subAteam', 'sky_tobs', 'sky_l', 'sky_m');
