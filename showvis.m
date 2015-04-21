% Script to show raw visibilities timeseries.
% pep/02Apr15
% Arguments:
%	fname  : visibility binary file name
%	colrng : caxis color range.
%	offset : The record offset from which to start showing images.
%	skip   : Skip some number of records.
%	nrecs  : Number of images to show. -1 => show all images.
%
% Returns:
%         void

function showvis (fname, offset, skip, nrecs, colrng)
	
	Nelem = 288;
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('showbinimages: fid < 0! Quitting.');
		return;
	end;

	if (nrecs < 0)
		[nrecs, tmin, tmax, dt] = getnrecs (fname);
	end;
	fprintf (2, '<-- Showing %d records.\n', nrecs);
	
	[acc, t_obs, freq] = readms2float (fid, offset, -1, Nelem);

	ind = 0
	while (ind < nrecs)
		imagesc (10*log10(abs (acc) - diag(diag(acc)))); colorbar;
		xlabel ('Ant'); ylabel ('Ant');
		title (sprintf ('%s, %.2f\n', datestr(mjdsec2datenum(t_obs)), freq));
        drawnow();
		for jind = 1:skip
			[acc, t_obs, freq] = readms2float (fid, -1, -1, Nelem);
		end;
        [acc, t_obs, freq] = readms2float (fid, -1, -1, Nelem);
		ind = ind + skip+1;
	end;	
