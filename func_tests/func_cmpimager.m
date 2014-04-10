% Script to compare FFT and DFT imagers on a standard dataset
% Comparison is done on the following:
%	- Difference between FFT and DFT images,
%	- Difference images from consecutive timeslices.
% pep/10Apr14

function func_cmpimager (fname)
	if (isempty(fname))
		fname = '/Users/peeyush/WORK/AARTFAAC/Reobs/20Nov13/r01/SB002_LBA_OUTER_8b2sbr01_1ch_1_convcal.bin';
	end;

	fid = fopen (fname, 'rb');
	offset = 20;
	winsize = 64;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500;						% size of gridded visibility matrix
    uvpad = 512;					% specifies if any padding needs to be added
	load ('../poslocal_outer.mat', 'poslocal');
	uloc = meshgrid (poslocal (:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal (:,2)) - meshgrid (poslocal(:,2)).';
	
	% Initialize datastructures
	acc = zeros (288, 288, winsize);
	tobs = zeros (1, winsize);
	fftmap = zeros (uvpad, uvpad, winsize);
	dftmap = fftmap;
	fftdifmap = zeros (uvpad, uvpad, winsize/2);
	dftdifmap = zeros (uvpad, uvpad, winsize/2);

	% Read in a couple of timeslices
	for ind = 1:winsize
		[acc(:,:,ind), tobs(ind), freq] = readms2float (fid, -1, -1, 288);
	end;

	fprintf (2, 'Generating fft images.\n');
	parm.type = 'wijnholds';
	parm.duv = duv;
	parm.Nuv = Nuv;
	parm.uvpad = uvpad;
	parm.lim = 0;
	parm.pa = [0 0 0 0];
	for ind = 1:winsize
		% vispad = genvisgrid (acc, u, v, parm, 1);
		tmpacc = acc (:,:,ind);
   		[radecmap, fftmap(:,:,ind), calvis, l, m] = ... 
		  fft_imager_sjw_radec (tmpacc(:), uloc(:), vloc(:), ... 
					duv, Nuv, uvpad, tobs(ind), freq, 0);
		if (mod (ind,2) == 0)
			fftdifmap (:,:,ind/2) = fftmap(:,:,ind) - fftmap(:,:,ind-1);
		end;
	end;

	fprintf (2, 'Generating DFT images.\n');
	for ind = 1:winsize
		tmpacc = acc (:,:,ind);
		tmpacc (isnan(tmpacc)) = 0;
		dftmap (:,:,ind) = acm2skyimage (tmpacc, poslocal(:,1), poslocal(:,2), freq, l, m);
		if (mod (ind,2) == 0)
			dftdifmap (:,:,ind/2) = dftmap(:,:,ind) - dftmap(:,:,ind-1);
		end;
	end;

	fprintf (2, 'Saving variables..\n');
	save ('func_cmpimager.mat', 'l', 'm', 'fftmap', 'dftmap', 'fftdifmap', 'dftdifmap');
	fclose (fid);
