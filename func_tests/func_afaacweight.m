% Script to evaluate the AARTFAAC PSF with different tapers, while keeping the
% density function identical. No GCF is applied, hence PSFs are DFT imaged.
% pep/21Apr14

%% LBA_OUTER
% load ('../poslocal_outer.mat', 'poslocal', 'posITRF');

% arrayconf = 'lba_outer';
% posfname  = 'poslocal_outer.mat';
arrayconf = 'lba_inner';
posfname  = 'poslocal_inner.mat';
freq = 60000000;
flagant = [];
deb=1;

fprintf (2, 'Generating no taper, natural weighting case for %s.\n', arrayconf);
% No tapering, natural weighting, DFT imaging.
fname = sprintf ('%s_uni_notap', arrayconf);
tparm.type = 'notaper';
tparm.pa(1) = -1;
tparm.pa(2) = -1;
tparm.pa(3) = -1;
tparm.pa(4) = -1;
wparm = []; gparm = [];
[l,m,psf,weight, intap, outtap] = genarraypsf (posfname, flagant, freq,...
	tparm, wparm, gparm, deb);


% Get figure handles
hand = get (0, 'Children');
print (hand(2), strcat (fname,'_psf.png'), '-dpng');
print (hand(1), strcat (fname,'_weights.png'), '-dpng');

% OUTER TAPER alone with different HPBW
% for out_sig = 10:50:300 %% NOTE: For 5m wavelength, 60 lambda = extent of array!
for out_sig = 10:10:60
	close all;
	fprintf (2, 'Generating Gaussian OUTERtap, sigma %d wavelengths.\n', ...
			 out_sig); 
	fprintf (2, 'HPBW (rad/deg) = %f/%f\n', 0.37/out_sig, ...
			(0.37/out_sig)*(180/pi));
	fname = sprintf ('%s_uni_outtap_%.0f', arrayconf, out_sig);
	tparm.type = 'gaussian';
	tparm.pa(1) = 0; % NO inner taper yet!
	tparm.pa(2) = 0;
	tparm.pa(3) = out_sig; % In wavelengths, not meters.
	tparm.pa(4) = out_sig;
	tparm.minlambda = 0;
	tparm.maxmeters = 350;
	wparm = []; gparm = [];
	[l,m,psf,weight, intap, outtap] = genarraypsf (posfname, flagant, freq,...
		tparm, wparm, gparm, deb);
	print (hand(2), strcat (fname,'_psf.png'), '-dpng');
	print (hand(1), strcat (fname,'_weights.png'), '-dpng');
end;


%% INNER TAPER alone with different HPBW
%% NOTE: TAPER FWHM specified in wavelengths!
for in_sig = 2:5:40
	close all;
	fprintf (2, 'Generating Gaussian INNERtap, sigma %d wavelengths.\n', ...
			 in_sig); 
	fname = sprintf ('%s_uni_intap_uvoff10m_%03.0f', arrayconf, in_sig);
	tparm.type = 'gaussian';
	tparm.pa(1) = in_sig;  
	tparm.pa(2) = in_sig; 
	tparm.pa(3) = -1;%NO outer taper yet!
	tparm.pa(4) = -1;
	tparm.minlambda = 10;
	tparm.maxmeters = 350;
	wparm = []; gparm = [];
	[l,m,psf,weight, intap, outtap] = genarraypsf (posfname, flagant, freq,...
		tparm, wparm, gparm, deb);
	hand = get (0, 'Children');
	print (hand(2), strcat (fname,'_psf.png'), '-dpng');
	print (hand(1), strcat (fname,'_weights.png'), '-dpng');
end;



%% Natural weighting with different cell sizes
fname = sprintf ('%s_uni_notap', arrayconf);
tparm.type = 'notaper';
tparm.pa(1) = -1;
tparm.pa(2) = -1;
tparm.pa(3) = -1;
tparm.pa(4) = -1;
gparm = [];

for cellrad = 10:5:50
	close all;
	wparm.type = 'uniform'; 
	wparm.cellrad = cellrad; % Meters.
	
	[l,m,psf,weight, intap, outtap] = genarraypsf (posfname, flagant, freq,...
		tparm, wparm, gparm, deb);
	fname = sprintf ('%s_uni_%dm_notap', arrayconf, cellrad);
	hand = get (0, 'Children');
	print (hand(2), strcat (fname,'_psf.png'), '-dpng');
	print (hand(1), strcat (fname,'_weights.png'), '-dpng');
end

%% Gridding: gaussian GCF with different parameters.
tparm.type = 'notaper';
tparm.pa(1) = -1;
tparm.pa(2) = -1;
tparm.pa(3) = -1;
tparm.pa(4) = -1;
wparm = [];
gparm.type = 'gaussian';
gparm.duv = 0.5;  % Now in wavelengths
gparm.Nuv = 241;
gparm.lim = 1.5; % In wavelength
gparm.uvpad = 256; 
gparm.fft = 1;
gparm.pa(1) = 0.5; gparm.pa(2) = gparm.pa(1); %% sigx/sigy of gaussian, wavelength units.

for cellrad = 0.5:1:3.5
	close all;
	gparm.pa(1) = cellrad; gparm.pa(2) = gparm.pa(1); %% sigx/sigy of gaussian, wavelength units.
	[l,m,psf,weight, intap, outtap] = genarraypsf (posfname, flagant, freq,...
		tparm, wparm, gparm, deb);
	fname = sprintf ('%s_uni_GCF_%.1f_notap', arrayconf, cellrad);
	hand = get (0, 'Children');
	print (hand(2), strcat (fname,'_psf.png'), '-dpng');
	print (hand(1), strcat (fname,'_weights.png'), '-dpng');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimized PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LBA_OUTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimized PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSF generated with optimized parameters based on above.
%% INNER and OUTER Taper, with parameters determined over the previous scan.
tparm.type = 'Gaussian';
tparm.minlambda = 10;
tparm.maxmeters = 350;
tparm.pa(1) = 17;  
tparm.pa(2) = 17; 
tparm.pa(3) = 60
tparm.pa(4) = 60;

wparm.type = 'uniform';
wparm.cellrad = 10;

gparm.type = 'Gaussian';
gparm.duv = 0.5;
gparm.Nuv = 500;
gparm.lim= 1.5;
gparm.uvpad= 512;
gparm.fft= 1;
gparm.pa= [0.5 0.5]
[l,m,psf,weight, intap, outtap, uvdist] = genarraypsf ('poslocal.mat', flagant, freq,...
	tparm, wparm, gparm, deb);

% Generate sample images from real data
fname = '/Users/peeyush/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB004_LBA_OUTER_SPREAD_1ch_1_convcal.bin';
flagant = [51, 238, 239, 273];
fid = fopen (fname, 'rb');
[acc, tobs, freqobs] = readms2float (fid, 1, -1, 288);
fclose (fid);

% DFT image
load ('poslocal_outer.mat', 'poslocal');
l = linspace (-1, 1, 512); m = l;
antsel = ones (1, length (poslocal(:,1)));
antsel (flagant) = 0;
visflag = ones (size (acc));
for ind = 1:length (flagant)
	visflag (:, flagant(ind)) = 0;
	visflag (flagant(ind), :) = 0;
end;
rem_ants = length(acc)-length(flagant);
acc_flag = reshape (acc(visflag == 1), [rem_ants rem_ants]);
skymap = acm2skyimage(acc_flag, poslocal(antsel==1, 1), poslocal(antsel==1, 2), freqobs, l, m);
figure;
imagesc (abs(skymap)); colorbar; caxis ([0 1200]);
title ('DFT image with natural taper, no weighting');

% acc_opt = acc(:) .* intap .* outtap .* 1./weight;
acc_opt = acc(:) .* (1./weight);
acc_opt = reshape (acc_opt (visflag==1), [rem_ants rem_ants]);
skymap_opt = acm2skyimage(acc_opt, poslocal(antsel==1, 1), poslocal(antsel==1, 2), freqobs, l, m);
figure;
imagesc (abs(skymap_opt)); colorbar; caxis ([0 1200]);
title ('DFT image with inner/outer taper, uniform weighting');

% FFT image with default pillbox
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
gparm.type = 'pillbox';
gparm.duv = 0.5; % NOTE: This is now in wavelength units.
gparm.Nuv = 500;
gparm.lim= 1.5; % 'Dont care' for pillbox
gparm.uvpad= 512;
gparm.fft= 1;
gparm.pa= [0.5 0.5] % 'Dont care' for pillbox
[radecmap, img.map, calvis, img.l, img.m] = ... 
	fft_imager_sjw_radec (acc_flag(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], tobs, freqobs, 0);

% FFT image with optimized PSF params
acc_opt = acc(:) .* intap .* outtap .* 1./weight;
acc_opt = acc_opt (visflag==1);
[radecmap, opt.map, calvis, opt.l, opt.m] = ... 
	fft_imager_sjw_radec (acc_opt(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], tobs, freqobs, 0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimized PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LBA_INNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimized PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSF generated with optimized parameters based on above.
%% INNER and OUTER Taper, with parameters determined over the previous scan.
tparm.type = 'Gaussian';
tparm.minlambda = 10;
tparm.maxmeters = 350;
tparm.pa(1) = 17;  
tparm.pa(2) = 17; 
tparm.pa(3) = 60
tparm.pa(4) = 60;

wparm.type = 'uniform';
wparm.cellrad = 10;

gparm.type = 'Gaussian';
gparm.duv = 0.5;
gparm.Nuv = 500;
gparm.lim= 1.5;
gparm.uvpad= 512;
gparm.fft= 1;
gparm.pa= [0.5 0.5]
[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq,...
	tparm, wparm, gparm, deb);

% Generate sample images from real data
fname_in = '/Users/peeyush/WORK/AARTFAAC/Reobs/11Jul12/LBA_INNER_BAND60/SB002_LBA_INNER_BAND60_1ch.bin'
flagant = [51, 238, 239, 273];
fid = fopen (fname, 'rb');
[acc, tobs, freqobs] = readms2float (fid, 1, -1, 288);
fclose (fid);

% DFT image
load ('poslocal_outer.mat', 'poslocal');
l = linspace (-1, 1, 512); m = l;
antsel = ones (1, length (poslocal(:,1)));
antsel (flagant) = 0;
visflag = ones (size (acc));
for ind = 1:length (flagant)
	visflag (:, flagant(ind)) = 0;
	visflag (flagant(ind), :) = 0;
end;
rem_ants = length(acc)-length(flagant);
acc_flag = reshape (acc(visflag == 1), [rem_ants rem_ants]);
skymap = acm2skyimage(acc_flag, poslocal(antsel==1, 1), poslocal(antsel==1, 2), freqobs, l, m);
figure;
imagesc (abs(skymap)); colorbar; caxis ([0 1200]);
title ('DFT image with natural taper, no weighting');

% acc_opt = acc(:) .* intap .* outtap .* 1./weight;
acc_opt = acc(:) .* (1./weight);
acc_opt = reshape (acc_opt (visflag==1), [rem_ants rem_ants]);
skymap_opt = acm2skyimage(acc_opt, poslocal(antsel==1, 1), poslocal(antsel==1, 2), freqobs, l, m);
figure;
imagesc (abs(skymap_opt)); colorbar; caxis ([0 1200]);
title ('DFT image with inner/outer taper, uniform weighting');

% FFT image with default pillbox
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
gparm.type = 'pillbox';
gparm.duv = 0.5; % NOTE: This is now in wavelength units.
gparm.Nuv = 500;
gparm.lim= 1.5; % 'Dont care' for pillbox
gparm.uvpad= 512;
gparm.fft= 1;
gparm.pa= [0.5 0.5] % 'Dont care' for pillbox
[radecmap, img.map, calvis, img.l, img.m] = ... 
	fft_imager_sjw_radec (acc_flag(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], tobs, freqobs, 0);

% FFT image with optimized PSF params
acc_opt = acc(:) .* intap .* outtap .* 1./weight;
acc_opt = acc_opt (visflag==1);
[radecmap, opt.map, calvis, opt.l, opt.m] = ... 
	fft_imager_sjw_radec (acc_opt(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], tobs, freqobs, 0);
