% Test driver script for generating gridded visibilities
% Test of genvisgrid.m
% pep/08Apr14

function func_testvisgrid ()
	% Open a calibrated vis. set
	fid = fopen ('~/WORK/AARTFAAC/Reobs/20Nov13/r01/SB002_LBA_OUTER_8b2sbr01_1ch_1_convcal.bin', 'r');
   [acc, tacc, freq] = readms2float (fid, 1, -1, 288);
	flagant = [129, 140, 149];

	% Create the u,v coordinates for the visibilities
	load ('../poslocal_outer.mat', 'poslocal');
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 

	% Generate the parameter structure for visgrid
	parm.type = 'Gaussian';
	parm.duv = 2.5;
	parm.Nuv = 240;
	parm.lim = 15; 
	parm.pa(1) = 5; parm.pa(2) = 5; % in meters


	% Generate gridded visibilities
	gridvis = genvisgrid (acc(:), uloc(:), vloc(:), parm, 1);
