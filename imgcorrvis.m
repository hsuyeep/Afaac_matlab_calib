% Script to generate stokes-I FITS images from GPU correlated visibilities available
% in a .vis file generated by the GPU correlator.

% Argument:
%	fname	: Array of filenames, one per subband, of binary .vis files.
%	obs		: Structure with observation related parameters. 
%		obs.npol	: Number of pols in input. [Default 2].
%		obs.nchan	: Number of chans in input.[Default 63].
%		obs.nelem	: Number of antenna elements in input.[Default 288].
%		obs.subband : Subband number of each file.
%		obs.freq	: Freq of each subband, in Hz.
%		obs.stokes	: Stokes required (I, Q, XX, YY) [Default I].
%		obs.bwidth	: Spectral integration needed on observation. Specify as a
%		  			  range of subbands [start:skip:end] [Default all presented
%		  			  subbands].
%		obs.imgspectint: Bool indicating whether spectral integration should occur
%			   		  pre or post imaging [Default postimaging].
%		obs.fov		: FoV to image, in deg[Default 180, all-sky].
%		obs.skip	: Stride of time records to image. Set as a python range
%					  [start:skip:end, -1 for end] [Default skip none].
%		obs.flagant : Preflagging of antennas;
%	fout	: Output folder, which will contain the generated FITS image files
%			  [Default same folder as that containing .vis files].

function [] = imgcorrvis (fname, obs, fout)

	% Validate the presence of the .vis file
	for i = 1:length (fname)
		assert (exist (fname{i}) == 2);
	end;
	nsub = length (fname);
    fprintf (2, '<-- Working with %d subbands.\n', nsub);
	
	% Setup the remainder of the obs structure, if uninitialized.
	if (isfield (obs, 'npol'  ) == 0) obs.npol  =   2; end;
	if (isfield (obs, 'nchan' ) == 0) obs.nchan =  63; end;
	if (isfield (obs, 'nelem' ) == 0) obs.nelem = 288; end;
	if (isfield (obs, 'nbline') == 0) obs.nbline= 41616; end;
	if (isfield (obs, 'stokes') == 0) obs.stokes=   4; end; % [XX=0,YY=1,XY=2,YX=3,I=4,Q=5]
	if (isfield (obs, 'fov'   ) == 0) obs.fov   = 180; end;
	if (isfield (obs, 'cal'   ) == 0) obs.cal   =   1; end;
	if (isfield (obs, 'skip'  ) == 0) obs.skip  =   1; end;
	if (isfield (obs, 'deb'   ) == 0) obs.deb   =   0; end;
	if (isfield (obs, 'ptSun' ) == 0) obs.ptSun =   1; end;
	if (isfield (obs, 'uvflag_x') == 0) obs.uvflag_x =   eye(obs.nelem); end;
	if (isfield (obs, 'uvflag_y') == 0) obs.uvflag_y =   eye(obs.nelem); end;
	if (isfield (obs, 'imgspectint') == 0) obs.imgspectint= 1; end;
	if (isfield (obs, 'bwidth') == 0) obs.bwidth= [1:nsub]; end;
	if (isfield (obs, 'flagant_x') == 0) obs.flagant_x = [18, 84, 142, 168, 262]; end; % For 8sb observations
	if (isfield (obs, 'flagant_y') == 0) obs.flagant_y = [18, 84, 142, 168, 262]; end;
	if (isfield (obs, 'posfilename') == 0) obs.posfilename = 'poslocal_outer.mat'; end;
	if (isfield (obs, 'sub') == 0) 
        for i = 1:length (fname) obs.sub(i) = 294+i; end;
    end;
	for i = 1:length(obs.sub)
			obs.freq(i) = (obs.sub(i)*195312.5);
	end;
	if (isfield (obs, 'gridparm') == 0)
		obs.gridparm.type = 'pillbox';
	    obs.gridparm.lim  = 0;

		% Control imaged field of view, independently of frequency.
	    % obs.gridparm.duv  = (0.5/180)*obs.fov;
	    obs.gridparm.duv  = 0.5;

	    obs.gridparm.Nuv  = 1024;	% size of gridded visibility matrix
	    obs.gridparm.uvpad= 1024;	% specifies if any padding needs to be added
	    obs.gridparm.fft  = 1;
	end;

    % Local horizon based coordinates of array in ITRF
    load (obs.posfilename, 'posITRF', 'poslocal'); 

    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

    [uloc_x, vloc_x] = gen_flagged_uvloc (uloc, vloc, obs.flagant_x);
    [uloc_y, vloc_y] = gen_flagged_uvloc (uloc, vloc, obs.flagant_y);

    info.nelem = obs.nelem;
    info.npol  = obs.npol;
    info.nchan = obs.nchan;
	for i = [1:nsub]
        info.freq = obs.freq(i);
        fprintf (1, '<-- Creating VisRec object for file %s.\n', fname{i});
		sbrecobj (i) = VisRec(fname{i}, info); 
	end;

	% Determine the number of records to be read, based on the skip parameter.
    % chan = [10, obs.nchan];
    chan = [10, 15];
	% Main loop handing the data.
    acm_x = zeros (nsub, obs.nelem, obs.nelem);
    acm_y = zeros (nsub, obs.nelem, obs.nelem);
    acm_tmp = zeros (obs.nelem);
    t1 = triu (ones (obs.nelem));
	for rec = 1:10
		for sb = 1:nsub
			% Read in a single record. Data stored internally.
			sbrecobj(sb).readRec ([1,0,0,1], [chan(1):chan(2)]);
		
			% Calibrate XX if stokes I or XX image is required
			if (obs.stokes >= 2 || obs.stokes == 0)
				% Average over selected channels, create nelem x nelem matrix.
                acm_tmp (t1 == 1) = mean (sbrecobj(sb).xx, 1);
                acm_tmp = acm_tmp + acm_tmp';
				acm_x(sb,:,:) = acm_tmp;

                if (obs.cal == 1)
	                fprintf (2, '\n<-- Calibrating XX for subband %d.\n', sb);
					sol_x(sb) = pelican_sunAteamsub (squeeze(acm_x(sb,:,:)), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_x, obs.flagant_x, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);
                    fprintf (1, '<-- Sigmas : %.5f\n', sol_x(sb).sigmas);
                end;
			end;

			% Calibrate YY if stokes I or YY image is required
			if (obs.stokes >= 2 || obs.stokes == 1)
                acm_tmp (t1 == 1) = mean (sbrecobj(sb).yy, 1);
                acm_tmp = acm_tmp + acm_tmp';
				acm_y(sb,:,:) = acm_tmp;
                
                if (obs.cal == 1)
                    fprintf (2, '\n<-- Calibrating YY for subband %d.\n', sb);
				    sol_y(sb) = pelican_sunAteamsub (squeeze(acm_y(sb,:,:)), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_y, obs.flagant_y, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);
                    fprintf (1, '<-- Sigmas : %.5f\n', sol_y(sb).sigmas);
                end;
			end;

		end;

		% If requested, splat all visibilities on a common grid
		if (obs.imgspectint == 0)

			% Splat all calibrated subbands on a common grid.
			integvis_x = zeros (size (sol_x.calvis));
			integvis_y = zeros (size (sol_y.calvis));
			% NOTE: This is probably incorrect! DONT USE NOW!
			for sb = 1:nsub
				[gridvis, gridviscnt, padgridrng, kern] = genvisgrid (sol_x.calvis, uloc_x, vloc_x, obs.gridparm, obs.freq(sb), 0);
				integvis_x = integvis_x + gridvis;
				[gridvis, gridviscnt, padgridrng, kern] = genvisgrid (sol_y.calvis, uloc_y, vloc_y, obs.gridparm, obs.freq(sb), 0);
				integvis_y = integvis_y + gridvis;
			end;

			% Image the commonly gridded visibilities. (Maybe just do a 2dFFT
			% right here)
   			[radecmap, integmap_x, calvis, l, m] = ... 
			    fft_imager_sjw_radec (integvis_x(:), uloc_x(:), vloc_x(:), ... 
					obs.gridparm, [], [], sbrecobj(1).trecstart, mean (obs.freq), 0);

   			[radecmap, integmap_y, calvis, l, m] = ... 
			    fft_imager_sjw_radec (integvis_y(:), uloc_y(:), vloc_y(:), ... 
					obs.gridparm, [], [], sbrecobj(1).trecstart, mean (obs.freq), 0);
		else
			% Image each subband separately
			for sb = 1:nsub

                fprintf (1, '\n<-- Imaging subband %d...\n', obs.sub(sb));
				if (obs.stokes >= 2 || obs.stokes == 0)
                    if (obs.cal == 1)
	   				    [radecmap, map_x(sb,:,:), calvis(sb,:,:), l, m] = ... 
					      fft_imager_sjw_radec (sol_x(sb).calvis(:), uloc_x(:), vloc_x(:), ... 
					    	obs.gridparm, [], [], sbrecobj(sb).trecstart, obs.freq(sb), 0);
                    else 
	   				    [radecmap, map_x(sb,:,:), calvis(sb,:,:), l, m] = ... 
					      fft_imager_sjw_radec (squeeze(acm_x(sb, :, :)), uloc(:), vloc(:), ... 
					    	obs.gridparm, [], [], sbrecobj(sb).trecstart, obs.freq(sb), 0);
                    end;
	            end;
	
				if (obs.stokes >= 2 || obs.stokes == 1)
                    if (obs.cal == 1)
	   				    [radecmap, map_y(sb,:,:), calvis(sb,:,:), l, m] = ... 
					         fft_imager_sjw_radec (sol_y(sb).calvis(:), uloc_y(:), vloc_y(:), ... 
						        obs.gridparm, [], [], sbrecobj(sb).trecstart, obs.freq(sb), 0);
                    else
	   				    [radecmap, map_y(sb,:,:), calvis(sb,:,:), l, m] = ... 
					         fft_imager_sjw_radec (squeeze(acm_y(sb,:,:)), uloc(:), vloc(:), ... 
						        obs.gridparm, [], [], sbrecobj(sb).trecstart, obs.freq(sb), 0);
                    end;
				end;
            end;
			
			% Add up the subband images to obtain the final integrated image.
			integmap = (squeeze(mean (map_x, 1) + mean (map_y, 1)))/2;
            % imagesc (squeeze(real(map_x(sb,:,:))));
            imagesc (squeeze(real(integmap)));
		end;

		% Generate FITS filename for this final image
		img.tobs 	  = sbrecobj(1).trecstart;
		img.pix2laxis = size (integmap, 1);
		img.pix2maxis = size (integmap, 2);
		img.freq 	  = mean (obs.freq);
		img.map 	  = integmap;
		
		integmapname = sprintf ('Sb%3d-%3d_R%02d-%02d_T%s.fits', obs.sub(1), obs.sub(end), chan(1), chan(2), datestr (mjdsec2datenum (sbrecobj(1).trecstart), 'dd-mm-yyyy_HH-MM-SS'));

		% Write out image as fits
		wrimg2fits (img, integmapname);
		% wrimg2bin(img, integmapname);
		fprintf (1, '<-- Writing out %s\n', integmapname);

	end % End of time dimension loop

end
