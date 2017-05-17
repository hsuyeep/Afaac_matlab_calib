% Function to read in a single binary visibility record, and return a
% subset of channels and pols. Also allows skipping time records.
% Arguments:
%	pol  : The pols to  be read in. Bool array, [XX, XY, YX, YY]	
%	chans: The channel subset to  be read in. Range [1:obj.nchans]	
% Returns:
%   flagdat : structure containg the actual data vectors for XX and YY pols, and
%		  metadata.
function obj = readRec(obj, pol, chans)
	assert (obj.fid > 0);

    stat = fseek (obj.fid, obj.recbytesize*obj.skip, 0);
    if (stat < 0)
        throw (MException ('VisRec.m:readRec()', ferror (obj.fid)));
    end;

	hdr = fread (obj.fid, 64, 'double'); % Read in the first 512 bytes
    if (feof (obj.fid) > 0)
        throw (MException ('VisRec.m:readRec()', ferror (obj.fid)));
    end;

	% Interpret the header of this data record. Hdr is first 512 bytes.
	meta = interpretHdr (obj, hdr);
    fprintf (1, '<-- Rec time: %s\n', datestr (mjdsec2datenum (meta.tstart)));
	
	assert (length (chans) <= obj.nchan);

	% Update internal state variables of this record
	obj.trecstart = meta.tstart;
	obj.trecend   = meta.tend;

	% dat = fread (obj.fid, obj.datfloatsize, [obj.nbline, obj.nchan, obj.npol, 2], 'single');
	dat = reshape (fread (obj.fid, obj.datfloatsize, 'single'), [2, obj.npol, obj.nchan, obj.nbline]);
    if (feof (obj.fid) > 0)
        throw (MException ('VisRec.m:readRec():', ferror (obj.fid)));
    end;

    if (obj.freqflag == 0)
        flagdat = mean (dat(:,:,chans,:), 3); % Carry out a blind averaging on the channel axis.
    else
        flagdat = obj.flagFreq (chans, dat);
    end;

	if (pol(1)) 
		obj.xx = complex (squeeze(flagdat(1,1,:)), squeeze(flagdat(2,1,:)));
	end;

	if (pol(4)) 
		obj.yy = complex (squeeze(flagdat(1,2,:)), squeeze(flagdat(2,2,:)));
	end;

    if (pol(2) || pol(3))
        fprintf (2, '### Ignoring XY and YX pols not present in incoming data. Ignoring.\n');
    end;

    dat2acm(obj);
endfunction;
