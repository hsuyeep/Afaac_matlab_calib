% Script to generate a subimage around selected sources via least squares
% imaging of every pixel. Uses calibrated visibilities.
% pep/17Jul13
% Arguments:
%	fname : File containing a timeseries of calibrated visibilities.
%	res   : Resolution of subimage in radians.
%  imsize : Size of subimage, assumed square, in pixels at above resolution.
%  srclist: indices in the 3CR srclist3CR of sources to image
% rodata  : Structure containing some required constants
%   debug : enables debug messages

function img = lsimg (fname, res, imsize, posfilename, debug)
		
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('func_calfluxrat: fid < 0! Quitting.');
		return;
	end;
		
	% Load various meta data
    rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
    rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
    rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	load (posfilename, 'posITRF', 'poslocal');
    rodata.posITRF_fl = posITRF; 	% NOTE: Not removing flagged vis. for now.
    rodata.poslocal = poslocal;
	load 'srclist3CR.mat';
	epoch = true;

	% NOTE: Current list assumes 3CR srclist3CR!
	% calsrcs = 3C[295, 219, 338, 410.1, 10, 84, 338]
	callist = [200, 133, 237]; %, 218, 4, 54, 237];
	ncal = length (callist);

	% Currently working only for a single timeslice.
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);

    tobs_jd = tobs/86400 + 2400000.5; % Convert from MJD secs. to JD day units

	fprintf (1, 'Time: %.2f, Freq: %.2f\n', tobs, freq);
	% Extract out calibrators within FoV and above horizon
	racal  = [srclist3CR(callist).alpha];
	deccal = [srclist3CR(callist).delta];
	fluxcal= [srclist3CR(callist).flux];
	srcposITRF = radectoITRF (racal, deccal, epoch, tobs_jd);
	upcals = srcposITRF * rodata.normal > 0;
	nsrc = sum (upcals);
	srcind = callist (upcals);

	% Spherical position of visible sources
	raupsrc = racal (upcals);	
	decupsrc = deccal (upcals); 

	% Filter out all flagged visbilities
	numnans = sum (sum (isnan(acc))); % List the number of usable baselines.
	acc (isnan(acc)) = 0; % Convert NaNs to 0s before using acc;

	img = zeros (nsrc, imsize+1, imsize+1);
	for src = 1:nsrc

		[lcen,mcen] = radectolm (raupsrc(src), decupsrc(src), tobs_jd,...
								 rodata.lon, rodata.lat, true);
		fprintf (1, '\nWorking on source %d.(%s, l,m=%.4f, %.4f)\n', ...
				src, srclist3CR(srcind(src)).name, lcen, mcen);
		% Generate a matrix of positions around the srclist3CR position at the 
		% specified resolution
		ravec = [(raupsrc (src)-(imsize/2)*res):res:...
				 (raupsrc (src)+(imsize/2)*res)];
		decvec= [(decupsrc(src)-(imsize/2)*res):res:...
				 (decupsrc(src)+(imsize/2)*res)];
		[ra, dec] = meshgrid (ravec, decvec);

		% Convert radec to lm coords.
		[l,m] = radectolm (ravec, decvec, tobs_jd, rodata.lon, ... 
							   rodata.lat, true);
		fprintf (1, 'Ravec: %s\n', num2str(ravec, '%.4f '));
		fprintf (1, 'l    : %s\n', num2str(l', '%.4f '));
		fprintf (1, 'decvec: %s\n', num2str(decvec, '%.4f '));
		fprintf (1, 'm    : %s\n', num2str(m', '%.4f '));
			
		% Convert position to cartesian ITRF positions.
		% Note that we can work only with a vector of positions, so image 
		% row by row.
		fprintf (1, '--> Row ');
		for row = 1:size (ra, 1)
			srcposITRFrow = radectoITRF (ra(row,:), dec(row,:), epoch, tobs_jd);
		
			% Generate Array steering vectors for each position in subimage
			A = exp(-(2 * pi * 1i * freq / 299792458) * ... 
					(rodata.posITRF_fl * srcposITRFrow.'));

			% Estimate powers in all pixels of subimage.
			img(src, row, :) = numnans*real(((abs(A' * A).^2) \ ...
										khatrirao(conj(A), A)')*acc (:));
			fprintf (1, ' %d', row);
		end;

		% Actually image the same specified region using fft imaging.
		if (debug > 2)
			figure;
			if (debug > 3) subplot (121); end;
			imagesc (l,m, squeeze (img (src, :,:)));

			if (debug > 3)
				% Carry out DFT imaging.
				fftmap (src, :, :) = ...
				acm2skyimage (acc, poslocal(:,1), poslocal(:,2), freq, l, m);

				% For carrying out FFT imaging.
    			% uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
		   		% vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
   				% [radecmap, fftmap (src,:,:), calvis, l, m] = ... 
				% 	  fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
				% 			duv, Nuv, uvpad, img.tobs, img.freq, radec);
				
				subplot (122);
				% Show the subimage.
				imagesc (l, m, squeeze (fftmap (src, :,:)));
				pause;
			end;
		end;
	end;

	if (debug > 2)
		for src = 1:nsrc
			imagesc (l,m, img (src, :,:));
			% title ('Source %s', 
		end;
	end;
