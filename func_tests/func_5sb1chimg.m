% Script to generate a wideband image using calibrated visibilities integrated
% to a single channel per subband. This is a short-cut from func_5sbimg.m, as 
% calibrated visiblities can be directly read from file.
% pep/10Jul13

addpath ../
load ('poslocal_outer.mat', 'posITRF', 'poslocal');

% func_5sb1chimg ()
fpath=  '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/';
fnames = {'SB000_LBA_OUTER_4min_1ch_1_convcal.bin', 
		  'SB001_LBA_OUTER_4min_1ch_1_convcal.bin', 
		  'SB002_LBA_OUTER_4min_1ch_1_convcal.bin', 
		  'SB003_LBA_OUTER_4min_1ch_1_convcal.bin', 
		  'SB004_LBA_OUTER_4min_1ch_1_convcal.bin'}; 


integtime = 60; % secs


% Imaging related
	mosaic = 0;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500; %1000                % size of gridded visibility matrix
    uvpad = 512; %1024              % specifies if any padding needs to be added
	nfacet = 3;
	facetsize = 256;
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

% Find time instances in each file
for fil = 1:length (fnames)
	[ntimes, tmin, tmax, dt] = getnrecs (strcat (fpath, fnames{fil}));
	fprintf (1, 'Found %d timesliecs in file %s\n', ntimes, fnames{fil});
	fid(fil) = fopen (strcat (fpath, fnames{fil}), 'r');
end;


% Datastructures
acc = zeros (288, 288, 5); % Temporary hold for a snapshot ACM for each subband.

% Holds the time integrated ACM from all subbands, where each subband ACM has been rephased to the center of the time window.
tsb_reph_integ_acc = zeros (288);  


% Holds the time integrated image from all subbands, where each subband ACM has 
% been rephased to the center of the time window.
if mosaic == 0
	sbintegimage = zeros (uvpad, uvpad, 5);
	timintegimage = zeros (uvpad);
else 
	sbintegimage = zeros (768, 768, 5);
	timintegimage = zeros (768);
end;

% For every time instance
for ts = 1:integtime

	% Read in all subbands
	for sb = 1:5
		[acc(:,:,sb), tobs(ts), freq(sb)] = readms2float (fid(sb), -1, -1, 288);
	end ;

	% ---------- Image plane stacking, each image at chosen freq. ------------ %
	for sb = 1:5
		% Rephase each sb's acm to central time
		reph = rephasetime (acc(:,:,sb), tobs(ts), tobs (1) + ntimes/2, ...
							freq(sb), posITRF);

		% Image this timeslice at sb frequency
	
		if (mosaic == 0)
	   		[radecmap, image_sb, calvis, l, m] = ... 
			  fft_imager_sjw_radec (reph(:), uloc(:), vloc(:), ... 
						duv, Nuv, uvpad, tobs(ts), freq(sb), 0);
		else
			% Image current timeslice. Generate a mosaic image.
	   		[image_sb, l, m] = ... 
			    genmosaic(tsb_reph_integ_acc, tobs (ts), freq(sb), nfacet, facetsize, ...
							uloc(:), vloc(:), posITRF, 0);
		end;

		% time integrate each subband rephased image separately
		sbintegimage (:,:,sb) = sbintegimage (:,:,sb) + image_sb;

	end;

	% ---------- Visibility averaging, without spectral gridding ------------ %
	integ_freq = zeros (288);
	% Combine subbands
	for sb = 1:5
		integ_freq = integ_freq + acc(:,:,sb);
	end;
	integ_freq = integ_freq/5;

	% Rephase to center of time window
	fprintf (1, 'Rec %d: Rephasing %.2f to time %.2f\n', ts, tobs (ts), tobs(1) + ntimes/2);
	reph = rephasetime (integ_freq, tobs(ts), tobs (1) + ntimes/2, freq(sb), ...
						posITRF);
	% integrate rephased acms in time.
	tsb_reph_integ_acc = tsb_reph_integ_acc + reph;
end;
	tsb_reph_integ_acc = tsb_reph_integ_acc / integtime;

	freq_integ = mean (freq);
	tim_integ = mean (tobs);

% Image the integrated acm
	lambda = 299792458/freq_integ; 		% in m.
    dl = (299792458/(freq_integ * uvpad * duv)); % dimensionless, in dircosunits
    lmax = dl * uvpad / 2;

	if (mosaic == 0)
		% Image current timeslice. Generate a zenith image.
   		[radecmap, integ_map, calvis, l, m] = ... 
		  fft_imager_sjw_radec (tsb_reph_integ_acc(:), uloc(:), vloc(:), ... 
					duv, Nuv, uvpad, tim_integ, freq_integ, 0);
	else
		% Image current timeslice. Generate a mosaic image.
   		[integ_map, l, m] = ... 
		    genmosaic(tsb_reph_integ_acc, tim_integ, freq_integ, nfacet, facetsize, ...
						uloc(:), vloc(:), posITRF, 0);
	end;
	mask = zeros (size (integ_map));
	mask (meshgrid (l).^2 + meshgrid(m).'.^2 < 1) = 1;

	imagesc (l, m, real (integ_map.*mask)); colorbar;
	title (sprintf ('%dsec, %d subband image. Freq. integ. ACM is rephased and time integ.', integtime, 5));
	xlabel ('l'); ylabel ('m');

	% Now show time averaged subband images, and stack them together
	sbintegimage = sbintegimage/integtime;
	figure;	
	for sb = 1:5
		% subplot (2,3,sb);
		% imagesc (real (sbintegimage (:,:,sb)));
		timintegimage = timintegimage + sbintegimage (:,:,sb);		
	end;
	timintegimage = timintegimage /5;
	% subplot (2,3,6);
	imagesc (l,m, real (timintegimage)); colorbar;
	title (sprintf ('%dsec, %d subband image. Freq. integ. ACM is rephased and time integ.', integtime, 5));
	xlabel ('l'); ylabel ('m');
