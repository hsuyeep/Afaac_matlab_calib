% Class to handle records from binary .vis files, as generated by the GPU correlator.
% pep/03Feb16

classdef VisRec < handle 
	% Class to abstract access to visibilities generated by the GPU correlator.
	% Supports the different polarizations, channels, antenna combinations
	% possible, represents a single timeslice of data.
	properties
		nchan       = 0;  % Number of channels in this record
		npol		= 0;  % Number of polarizations
		nelem	    = 0;  % Number of antenna elements
		nbline      = 0;  % Number of baselines
		freq		= 0;  % Frequency of the subband of this visibility set.
		fname		= {}; % Filename of file containing data
		fid			= 0;  % fid of file.
		currec		= 0;  % Current record number in the binary file.
		tfilestart  = 0;  % Timestamp of first record in binary file.
		tfileend	= 0;  % Timestamp of last record in binary file.
		trecstart   = 0;  % Timestamp of first record in binary file.
		trecend	    = 0;  % Timestamp of last record in binary file.
		skip		= 0;  % Number of records to skip while reading data.
		deb			= 0;  % Debug level, for verbosity of messages.
        recbytesize = 0;
        datfloatsize= 0;
        xx          = [];
        xy          = [];
        yx          = [];
        yy          = [];

	end; % End of properties.

	methods 

		% Constructor
		% Opens the binary file, initializes internal variables.
		% If metainfo. cannot be determined from the data stream,
		% uses the externally supplied values.
        % Arguments:
        %   fname   : Filename of visibility binary file
        %   info    : Structure containing observational information.
		function  obj = VisRec (fname, info)
			assert (exist (fname) == 2);
			obj.fname = fname;
			% try
			obj.fid  = fopen (obj.fname, 'rb');	
			% catch
		    % 	fprintf (2, '### Error in opening file %s.\n', obj.fname);
			% end;
						
			% obj.recdat = RecDat(); % Self contained data structure to hold a single records' data.
			obj.getVisMeta (info); % Initialize internal data based on input data.
		end;

		% Destructor
		function delete (obj)
			if (obj.fid > 0)
				fclose (obj.fid);
			end;
		end;

		% Examine data stream, generate metadata structure.
		% Called only once to initialize internal data.
		% Arguments:
		%	info	: Structure containing the preinitalized parameters of this
		%			  observation. Note: raw visibility header contains only the
		%			  magic, start and end times, nothing more!
		function info = getVisMeta (obj, info)
			fprintf (1, '<-- Generating observation scan visibility meta data...\n');
			assert (obj.fid > 0);
			hdr = fread (obj.fid, 64, 'double');
			meta = obj.interpretHdr (hdr);
			if (meta.magic == 999878658) % = 0x3B98F002, correlator magic number.
				fprintf (1, '<-- Found Magic number 0x%x in first record.\n',meta.magic);		
				obj.tfilestart = meta.tstart;
				fprintf (1, '<-- First record start time: %s.\n', datestr (mjdsec2datenum(obj.tfilestart)));
				obj.nchan = info.nchan;
				obj.nelem = info.nelem;
				obj.npol  = info.npol;
				obj.freq  = info.freq;
                obj.nbline= info.nelem*(info.nelem+1)/2;
				% In units of floats, the native output datatype of the
				% correlator.
				obj.datfloatsize = 2*obj.npol*obj.nchan*obj.nelem*(obj.nelem+1)/2; 
				obj.recbytesize = 512 + obj.datfloatsize*4; % Note: In bytes
			end;
		
            fprintf (1, '<-- Record size : %d Bytes, datasize: %d floats.\n', obj.recbytesize, obj.datfloatsize);

			% Move to end of the file to obtain last record.
%			fseek (obj.fid, -(obj.recbytesize+1), 'eof');
%			hdr = fread (obj.fid, 64, 'double'); % Read in hdr of the last record.
%			meta = obj.interpretHdr (hdr);
%			if (meta.magic == 999878658)
%				fprintf (1, '<-- Found Magic number 0x%x in last record. %x\n', meta.magic);		
%				obj.tfileend = meta.tstart;
%				fprintf (1, '<-- Last record start time: %s.\n', datestr(mjdsec2datenum(obj.tfileend)));		
%			end;
			fseek (obj.fid, 0, 'bof'); % Move to the beginning of the file for further operations.
		end;
		
		% Function interprets the binary header and returns meta information.
		% Arguments:
		%	hdr : array of 512 uint8, corresponding to the first 512 bytes of
		%		  the visibility record.
		% Returns :
		%	meta: A structure containing metadata as extracted from the header.
		function meta = interpretHdr (obj, hdr)
			% Set the object magic value
            tmp = typecast (hdr(1), 'uint32'); 
			meta.magic  = tmp(1);
            
			meta.tstart = unixtime2mjdsec (hdr(2));
			meta.tend   = unixtime2mjdsec (hdr(3));
		end;

		% Arguments:
		%	pol  : The pols to  be read in. Bool array, [XX, XY, YX, YY]	
		%	chans: The channel subset to  be read in. Range [1:obj.nchans]	
		% Returns:
		%	dat : structure containg the actual data vectors for XX and YY pols, and
		%		  metadata.
		function dat = readRec(obj, pol, chans)
			assert (obj.fid > 0);
			hdr = fread (obj.fid, 64, 'double'); % Read in the first 512 bytes

			% Interpret the header of this data record. Hdr is first 512 bytes.
			meta = obj.interpretHdr (hdr);
			
			assert (length (chans) <= obj.nchan);

			% Update internal state variables of this record
%			obj.recdat.nchan  = obj.nchan;
%			obj.recdat.npol   = obj.npol;
%			obj.recdat.nelem  = obj.nelem;
%			obj.recdat.freq   = obj.freq;
			obj.trecstart = meta.tstart;
			obj.trecend   = meta.tend;
% 			obj.recdat.dt = meta.tend - meta.tstart;

			% dat = fread (obj.fid, obj.datfloatsize, [obj.nbline, obj.nchan, obj.npol, 2], 'single');
			dat = reshape (fread (obj.fid, obj.datfloatsize, 'single'), [2, obj.npol, obj.nchan, obj.nbline]);
			% Fill in the RecDat object with desired channels and pols.
			% The fastest index in the linear float data array is nbline, nchan,
			% npol, re/im
			if (pol(1)) 
				obj.xx = complex (squeeze(dat(1,1,chans,:)), squeeze(dat(2,1,chans,:)));
			end;
			if (pol(4)) 
				obj.yy = complex (squeeze(dat(1,2,chans,:)), squeeze(dat(2,2,chans,:)));
			end;
			% NOTE: Ignoring XY and YX pols for now.
		end;
			
	end % End of methods

end % End of classdef
