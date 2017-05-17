% Examine data stream, generate metadata structure.
% Called only once to initialize internal data.
% Arguments:
%	info	: Structure containing the preinitalized parameters of this
%			  observation. Note: raw visibility header contains only the
%			  magic, start and end times, nothing more!
function obj = getVisMeta (obj, info)
	fprintf (1, '<-- Generating observation scan visibility meta data...\n');
	assert (obj.fid > 0);
	hdr = fread (obj.fid, 64, 'double');
	meta = interpretHdr (obj, hdr);
    fprintf (1, '%d', meta.magic);
	if (meta.magic == 999878658) % = 0x3B98F002, correlator magic number.
		fprintf (1, '<-- Found Magic number 0x%x in first record.\n',meta.magic);		
        fprintf (1, '<-- File contains RAW Correlations.\n');
		obj.tfilestart = meta.tstart;
        obj.trecstart = meta.tstart;
		fprintf (1, '<-- First record start time: %s.\n', datestr (mjdsec2datenum(obj.tfilestart)));
        if (isfield (info, 'nchan')) obj.nchan = info.nchan; end;
		if (isfield (info, 'nelem')) obj.nelem = info.nelem; end;
		if (isfield (info, 'npol')) obj.npol  = info.npol; end;
		if (isfield (info, 'freq')) obj.freq  = info.freq; end;
		if (isfield (info, 'sub')) obj.freq  = info.sub*195312.5; end;
        if (isfield (info, 'nbline')) obj.nbline= info.nelem*(info.nelem+1)/2; end;
		% In units of floats, the native output datatype of the
		% correlator.
		obj.datfloatsize = 2*obj.npol*obj.nchan*obj.nbline;
		obj.recbytesize = 512 + obj.datfloatsize*4; % Note: In bytes

    % 0x4141525446414143, calibrated visibilities magic number.
    elseif (meta.magic == 0x41415254) % Looking at first 32 bits only
    % elseif (meta.magic == 0x4141525446414143)
        fprintf (1, '<-- Found Magic number 0x%x in first record.\n',meta.magic);
        fprintf (1, '<-- File contains CALIBRATED Correlations.\n');
        obj.tfilestart = meta.tstart;
        obj.trecstart = meta.tstart;
        fprintf (1, '<-- First record start time: %s.\n', datestr(datenum(obj.tfilestart)));
        if (isfield (info, 'nchan')) obj.nchan = info.nchan; end;
        if (isfield (info, 'nelem')) obj.nelem = info.nelem; end;
        if (isfield (info, 'npol')) obj.npol  = info.npol; end;
        if (isfield (info, 'freq')) obj.freq  = info.freq; end;
        if (isfield (info, 'sub')) obj.freq  = info.sub*195312.5; end;
        if (isfield (info, 'nbline'))
            obj.nbline= info.nelem*info.nelem;
        end;
        % In units of floats, the native output datatype of the
        % correlator.
        obj.datfloatsize = 2*obj.npol*obj.nchan*obj.nbline;
        obj.recbytesize = 512 + obj.datfloatsize*4; % Note: In bytes
	end;

    fprintf (1, '<-- Record size : %d Bytes, datasize: %d floats.\n', obj.recbytesize, obj.datfloatsize);

    % Move to next record to obtain integration time
	fseek (obj.fid, 0, 'bof'); % Move to the beginning of the file for further operations.
    fseek (obj.fid, obj.recbytesize, 'cof'); % Move past first record, which can have a 0 in its timestamp.
    
	hdr = fread (obj.fid, obj.recbytesize/8, 'double'); % Read in the full record as doubles.
	meta = interpretHdr (obj, hdr);
    obj.tfilestart = meta.tstart; % Rejecting first record.
    obj.trecstart = meta.tstart;
    prevtime = obj.tfilestart;
    for i = 1:3 % Arbitrarily reading 3 records to determine dt.
		hdr = fread (obj.fid, obj.recbytesize/8, 'double'); % Read in the full record as doubles.
		meta = interpretHdr (obj, hdr);
        obj.trecstart = meta.tstart;
        dt(i) = meta.tstart - prevtime;
        prevtime = meta.tstart;
		if (meta.magic == 999878658) % = 0x3B98F002, correlator magic number.
			fprintf (1, '<-- Found Magic number 0x%x in record %d, time %s, dt %fsec.\n', ...
                     meta.magic, i, datestr (mjdsec2datenum(obj.trecstart)), dt(i));		
        end;
    end;
    obj.dt = mean (dt);
        
    
	% Move to end of the file to obtain last record.
	fseek (obj.fid, -(obj.recbytesize), 'eof');
	hdr = fread (obj.fid, obj.recbytesize/8, 'double'); % Read in hdr of the last record.
	meta = interpretHdr (obj, hdr);
	if (meta.magic == 999878658)
		fprintf (1, '<-- Found Magic number 0x%x in last record. %x\n', meta.magic);		
		obj.tfileend = meta.tstart;
		fprintf (1, '<-- Last record start time: %s.\n', datestr(mjdsec2datenum(obj.tfileend)));		
	end;
	fseek (obj.fid, 0, 'bof'); % Move to the beginning of the file for further operations.
end;
