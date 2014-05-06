% Test driver script for generating gridded visibilities
% Test of genvisgrid.m
% pep/08Apr14

function func_testvisgrid ()

	% If a real image needs to be created.
	% Open a calibrated vis. set
	fid = fopen ('~/WORK/AARTFAAC/Reobs/20Nov13/r01/SB002_LBA_OUTER_8b2sbr01_1ch_1_convcal.bin', 'r');
   [acc, tacc, freq] = readms2float (fid, 1, -1, 288);
	flagant = [129, 140, 149];

	sampfn = ones (288); flagant = [];
	freq = 60000000; % LBA Resonance

	% Create the u,v coordinates for the visibilities
	load ('../poslocal_outer.mat', 'poslocal');
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 

	% Generate the parameter structure for visgrid
	parm.type = 'pillbox';
	% parm.type = 'Gaussian';
	parm.duv   = 0.5; % In wavelengths
	parm.Nuv   = 1001; % Size of gridded visibilities, no. of rows/cols.
	parm.uvpad = 1024;
	parm.lim   = 1.5; % In wavelength
	parm.pa(1) = 0.5; % 2D Gaussian x/ysigma, in wavelengths
	parm.pa(2) = 0.5; 
	parm.fft   =   0; % Carry out FFT imaging.   


	% Generate gridded visibilities with different 
	gridvis = genvisgrid (sampfn(:), uloc_flag(:), vloc_flag(:), parm, freq, 1);

	% Generate an image with GCF.
	% [rdsky, map, vispad, l, m] = fft_imager_sjw_radec (acc, uloc_flag, ...
% 						vloc_flag, parm, [], [], 0, freq, 0);
	
% 	imagesc (l, m, abs(map));
