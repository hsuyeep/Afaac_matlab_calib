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

%% PSF generated with optimized parameters based on above.
%% INNER and OUTER Taper, with parameters determined over the previous scan.
%	tparm.pa(1) = 60;  
%	tparm.pa(2) = 60; 
%	tparm.pa(3) = 7;%NO outer taper yet!
%	tparm.pa(4) = 17;
%	tparm.minlambda = 10;
%	tparm.maxmeters = 350;
%	wparm = []; gparm = [];
%	[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq,...
%		tparm, wparm, gparm, deb);
