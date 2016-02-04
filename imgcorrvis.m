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
		assert (exist (fname{i} == 0);
	end;
	nsubs = length (fname);
	
	% Setup the remainder of the obs structure, if uninitialized.
	if (isempty (obs.npol )) obs.npol  =   2; end;
	if (isempty (obs.nchan)) obs.nchan =  63; end;
	if (isempty (obs.nelem)) obs.nelem = 288; end;
	if (isempty (obs.stokes))obs.stokes=   0; end;
	if (isempty (obs.fov  )) obs.fov   = 180; end;
	if (isempty (obs.skip )) obs.skip  =   1; end;
	if (isempty (obs.deb  )) obs.deb   =   0; end;
	if (isempty (obs.ptSun)) obs.ptSun =   1; end;
	if (isempty (obs.uvflag_x)) obs.uvflag_x =   eye(obs.nelem); end;
	if (isempty (obs.uvflag_y)) obs.uvflag_y =   eye(obs.nelem); end;
	if (isempty (obs.imgspectint)) obs.imgspectint= 1; end;
	if (isempty (obs.bwidth)) obs.bwidth= [1:nsubs]; end;
	if (isempty (obs.flagant_x)) obs.flagant_x = []; end;
	if (isempty (obs.flagant_y)) obs.flagant_y = []; end;
	if (isempty (obs.posfile)) obs.posfile = 'poslocal_outer.mat'; end;
	if (isempty (obs.freq)) 
		for i = 1:length (fname)
			obs.freq(i) = (295+i)*195312.5;
		end;
	end;	
	if (isempty (obs.gridparm))
		obs.gridparm.type = 'pillbox';
	    obs.gridparm.lim  = 0;
	    obs.gridparm.duv  = 0.5;	% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
	    obs.gridparm.Nuv  = 1500;	% size of gridded visibility matrix
	    obs.gridparm.uvpad= 1536;	% specifies if any padding needs to be added
	    obs.gridparm.fft  = 1;
	end;

	% Imaging related setup
	lambda = 299792458/obs.freq;	% in m.
    dl = (299792458/(obs.freq * obs.gridparm.uvpad * obs.gridparm.duv)); % dimensionless, in dir. cos. units
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * obs.gridparm.uvpad / 2;

    % Local horizon based coordinates of array in ITRF
    load (posfilename, 'posITRF', 'poslocal'); 

    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

    [uloc_x, vloc_x] = gen_flagged_uvloc (uloc, vloc, obs.flagant_x);
    [uloc_y, vloc_y] = gen_flagged_uvloc (uloc, vloc, obs.flagant_y);

	% Validate the specified observational parameters with those specified in
	% the obs. structure.
	for i = [1:nsubs]
		sbrecobj (i) = VisRec(fname{i}, obs.freq(i));
	end;

	% Determine the number of records to be read, based on the skip parameter.

	% Main loop handing the data.
	for rec = 1:nrec
		for sb = 1:nsubs
			% Read in a single record. Data stored internally.
			sbrecobj(sb).readRec ([1,0,0,1], [1:obs.nchan]);
		
			% Calibrate XX and YY separately for each subband.
			if (obs.stokes >= 2 || obs.stokes == 0)
				% Average over selected channels, create nelem x nelem matrix.
				acm_x(sb) = mean (sbrecobj(sb).recdat.xx, 2);
				sol_x(sb) = pelican_sunAteamsub (acm_x(sb), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_x, obs.flagant_x, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);
			end;

			if (obs.stokes >= 2 || obs.stokes == 1)
				acm_y(sb) = mean (sbrecobj(sb).recdat.yy, 2);
				sol_y(sb) = pelican_sunAteamsub (acm_y(sb), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_y, obs.flagant_y, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);
			end;

		end;

		% If requested, splat all visibilities on a common grid
		if (imgspectint == 0)
			% Splat all calibrated subbands on a common grid.
			% TODO;

			% Image the commonly gridded visibilities.
   			[radecmap, integmap, calvis, l, m] = ... 
			  fft_imager_sjw_radec (sol_x.calvis(:), uloc_x(:), vloc_x(:), ... 
					obs.gridparm, [], [], sbrecobj(1).trecstart, sbrecobj(1).freq, 0);
		else
			% Image each subband separately
			for sb = 1:nsubs
   				[radecmap, map(sb), calvis(sb), l, m] = ... 
			  fft_imager_sjw_radec (sol_x(sb).calvis(:), uloc_x(:), vloc_x(:), ... 
					obs.gridparm, [], [], sbrecobj(sb).trecstart, sbrecobj(sb).freq, 0);
			end;
			
			% Add up the subband images to obtain the final integrated image.
			integmap = mean (map, 1);
		end;

		% Generate FITS filename for this final image
		integmapname = sprintf ('Int%2d_R%02d-%02d_T%s.fits', nsubs, chan(1), chan(2), datestr (mjdsec2datenum (sbrecobj(1).trecstart), 'dd-mm-yyyy_HH-MM-SS'));

		% Write out image as fits
		wrimg2fits (integmap, integmapname);
	end; % End of time dimension loop

end;
