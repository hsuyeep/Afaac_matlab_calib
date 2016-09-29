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
        freqflag    = 0;  % Bool to indicate whether channel axis flagging should be carried out.
        timeflag    = 0;  % Bool to indicate whether temporal flagging should be carried out.
		fname		= {}; % Filename of file containing data
		fid			= 0;  % fid of file.
		currec		= 0;  % Current record number in the binary file.
		tfilestart  = 0;  % Timestamp of first record in binary file.
		tfileend	= 0;  % Timestamp of last record in binary file.
		trecstart   = 0;  % Timestamp of start of currently read record.
		trecend	    = 0;  % Timestamp of end of currently read record.
		skip		= 0;  % Number of records to skip while reading data.
		deb			= 0;  % Debug level, for verbosity of messages.
        recbytesize = 0;
        datfloatsize= 0;
        flagsig     = 2.5;% Sigma threshold at which to clip
        xx          = [];
        xy          = [];
        yx          = [];
        yy          = [];
        acm_xx      = [];
        acm_yy      = [];
        dt          = 0;

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
			try
			    obj.fid  = fopen (obj.fname, 'rb');	
			catch ME
		     	fprintf (2, '### Error in opening file %s.\n', obj.fname);
			end;

			if (isfield (info, 'nchan' ) == 0) obj.nchan =  63; else obj.nchan = info.nchan; end;
			if (isfield (info, 'npol'  ) == 0) obj.npol  =   2; else obj.npol = info.npol; end;
			if (isfield (info, 'nelem' ) == 0) obj.nelem = 288; else obj.nelem = info.nelem; end;
			if (isfield (info, 'nbline') == 0) obj.nbline= 41616; else obj.nbline = info.nbline; end;
			if (isfield (info, 'freqflag') == 0) obj.freqflag= 0; else obj.freqflag = info.freqflag; end;
			if (isfield (info, 'timeflag') == 0) obj.timeflag= 0; else obj.timeflag = info.timeflag; end;
			if (isfield (info, 'freq') == 0) obj.freq = 0; else obj.freq = info.freq; end;
			if (isfield (info, 'skip'  ) == 0) obj.skip  =   0; else obj.skip = info.skip; end;
			if (isfield (info, 'deb'   ) == 0) obj.deb   =   0; else obj.deb = info.deb; end;
						
			% obj.recdat = RecDat(); % Self contained data structure to hold a single records' data.
			obj.getVisMeta (info); % Initialize internal data based on input data.
		end;


		% Destructor
		function delete (obj)
			if (obj.fid > 0)
				fclose (obj.fid);
			end;
		end;


        % Return the internal state of the VisRec object.
        function [obs] = getObsInfo (obj)
		    if (isfield (obj, 'nchan' ) == 0)    obs.nchan     = obj.nchan; end;    
		    if (isfield (obj, 'npol'  ) == 0)    obs.npol      = obj.npol; end;     
		    if (isfield (obj, 'nelem' ) == 0)    obs.nelem     = obj.nelem; end;    
		    if (isfield (obj, 'nbline') == 0)    obs.nbline    = obj.nbline; end;   
		    if (isfield (obj, 'freqflag')== 0)   obs.freqflag  = obj.freqflag; end; 
		    if (isfield (obj, 'timeflag')== 0)   obs.timeflag  = obj.timeflag; end; 
		    if (isfield (obj, 'freq')   == 0)    obs.freq      = obj.freq; end; 
		    if (isfield (obj, 'skip')   == 0)    obs.skip      = obj.skip; end; 
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
				obj.datfloatsize = 2*obj.npol*obj.nchan*obj.nelem*(obj.nelem+1)/2; 
				obj.recbytesize = 512 + obj.datfloatsize*4; % Note: In bytes
			end;
		
            fprintf (1, '<-- Record size : %d Bytes, datasize: %d floats.\n', obj.recbytesize, obj.datfloatsize);

            % Move to next record to obtain integration time
			fseek (obj.fid, 0, 'bof'); % Move to the beginning of the file for further operations.
            fseek (obj.fid, obj.recbytesize, 'cof'); % Move past first record, which can have a 0 in its timestamp.
            
    		hdr = fread (obj.fid, obj.recbytesize/8, 'double'); % Read in the full record as doubles.
    		meta = obj.interpretHdr (hdr);
            obj.tfilestart = meta.tstart; % Rejecting first record.
            obj.trecstart = meta.tstart;
            prevtime = obj.tfilestart;
            for i = 1:3 % Arbitrarily reading 3 records to determine dt.
    			hdr = fread (obj.fid, obj.recbytesize/8, 'double'); % Read in the full record as doubles.
    			meta = obj.interpretHdr (hdr);
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
			meta = obj.interpretHdr (hdr);
			if (meta.magic == 999878658)
				fprintf (1, '<-- Found Magic number 0x%x in last record. %x\n', meta.magic);		
				obj.tfileend = meta.tstart;
				fprintf (1, '<-- Last record start time: %s.\n', datestr(mjdsec2datenum(obj.tfileend)));		
			end;
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


        % Function to carry out visibility splatting on a common grid, but
        % across frequencies.
        % Arguments:
        %   vis     : Visibilities over different subbands/pols.
        %  freq     : Frequencies of the various subbands.
        %   gparm   : Gridding parameter structure.
        % poslocal  : Locational information on the antennas.
        % Returns:
        %   gvis    : Gridded visibilities.
        function vis = gridFreq (obj, vis, freq, gparm)

        end;


        % Function to carry out visibility based flagging along the freq.
        % axis. Statistics are generated over the channel axis.
        % Arguments:
        %   chans   : The Channels over which to calculate statistics.
        %   dat     : Reshaped float array of visibilities.
        % Returns:
        %   flagdat : visibilities mean after ignoring flagged channels.
        function flagdat = flagFreq (obj, chans, dat)
            flagdat = zeros (2, obj.npol, obj.nbline);

            % Doing it the inefficient way for now.
            for v = 1:obj.nbline
                for p = 1:obj.npol
                    for r = 1:2 % re/im
                        seldat = squeeze (dat (r, p, chans, v));
                        sel = ones (obj.nchan, 1) & ~isnan (seldat);
                        for i = 1:5
			                mvis = mean (seldat(sel == 1));
			                svis = std  (seldat(sel == 1));
			                sel  = (abs (seldat(sel) - mvis) < obj.flagsig*svis);
                        end
                        flagdat (r,p,v) = mean (seldat(sel));
                    end;
                end;
            end;
            
        end;


        % Function to convert data stored in the xx/yy linear arrays into
        % a square symmetric Array Covariance Matrix.
        function dat2acm (obj)
            assert (~isempty (obj.xx));
            assert (~isempty (obj.yy));
            t = triu (ones(obj.nelem));
            obj.acm_xx = zeros (obj.nelem);
            obj.acm_xx (t == 1) = obj.xx; 
            acm_diag = diag (diag (obj.acm_xx)); 
            obj.acm_xx = conj (obj.acm_xx + obj.acm_xx' - acm_diag);

            obj.acm_yy = zeros (obj.nelem);
            obj.acm_yy (t == 1) = obj.yy; 
            acm_diag = diag (diag (obj.acm_yy)); 
            obj.acm_yy = conj (obj.acm_yy + obj.acm_yy' - acm_diag);
        end;

        % Function to skip raw records in the .vis file.
        % Assumes that the observation parameters of data are already
        % initialized.
        % Arguments:
        %   rec : Number of records to skip.
        %   unit: string ('rec', 'dt', 'mjdsec') Whether to skip records or seconds, or to
        %   an absolute time.
        % Returns:
        %   Void.
        function skipRec (obj, rec, unit)
            % If unit is time (sec), determine number of records to skip based
            % on integration time and current time.
            switch lower(unit)
                case 'sec' % Should be in seconds
                    recs2move = int32(rec/obj.dt);
                    fprintf (1, '<-- Moving %d recs for %d seconds.\n', recs2move, rec)
                    stat = fseek (obj.fid, obj.recbytesize*recs2move, 0);
                    if (stat < 0)
                        throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
                    end;
                
                case 'mjdsec'
                    assert (obj.tfileend > rec);
                    assert (obj.tfilestart < rec);
                    recs2move = int32((rec-obj.trecstart)/obj.dt);
                    fprintf (1, '<-- Moving %d recs (class %s) for %f seconds.\n', recs2move, class (recs2move), rec - obj.trecstart)
                    stat = fseek (obj.fid, obj.recbytesize*recs2move, 0);
                    if (stat < 0)
                        throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
                    end;
                
                case 'rec' 
                    stat = fseek (obj.fid, obj.recbytesize*rec, 0);
                    if (stat < 0)
                        throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
                    end;
            
                otherwise
                    fprintf (2, '### Unknown skip option %s.\n', unit);
            end;
        end;

        % Function to read in a single binary visibility record, and return a
        % subset of channels and pols. Also allows skipping time records.
		% Arguments:
		%	pol  : The pols to  be read in. Bool array, [XX, XY, YX, YY]	
		%	chans: The channel subset to  be read in. Range [1:obj.nchans]	
		% Returns:
		%   flagdat : structure containg the actual data vectors for XX and YY pols, and
		%		  metadata.
		function dat = readRec(obj, pol, chans)
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
			meta = obj.interpretHdr (hdr);
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
			% NOTE: Ignoring XY and YX pols for now.
		end;


        % Function to read in multiple timeslices of data, sigmaclip in time
        % axis and present an averaged visibility to the user. For implementing
        % time domain averaging in visibilities before imaging.
        % Calls readRec() underneath.
		% Arguments:
		%	pol  : The pols to  be read in. Bool array, [XX, XY, YX, YY]	
		%	chans: The channel subset to  be read in. Range [1:obj.nchans]	
        % tslices: The number of timeslices to average over.
		% Returns:
		%   dat  : Time averaged visibilities, flagged both along the channel and time
		%          axis.
		function dat = readTimeAvgRec (obj, pol, chans, tslices)
            
            for ind = 1:tslices
                obj.readRec (pol, chans);
                xxtslices(ind, :) = obj.xx;
                yytslices(ind, :) = obj.yy;
            end;

            if (obj.timeflag)
	            % Selection mask
	            selxx = ones (size (xxtslices));
	            selyy = ones (size (yytslices));
                % ### UNIMPLEMENTED, DONT USE!
	            while 1
	                m = nanmean (xxtslices);
	                s = nanstd  (xxtslices);
	                xxtslices(abs(xxtslices-repmat (m, [obj.nbline, 1])) > repmat (obj.flagsig*abs(s), [obj.nbline, 1])) = nan;
	
	                m = nanmean (yytslices);
	                s = nanstd  (yytslices);
	                yytslices(abs(yytslices-repmat (m, [obj.nbline, 1])) > repmat (obj.flagsig*abs(s), [obj.nbline, 1])) = nan;
	            end;
            else
                dat.xx = mean (xxtslices, 1);
                dat.yy = mean (yytslices, 1);
            end;
        end;


        % Function to generate a calibrated map from the observed
        % visibilities. Returns generated images in the img structure
        % Both XX and YY images are returned.
        %
        % Arguments:
        %   arrayconfig : One of 'lba_outer', 'lba_inner'
        %   flagant     : Index of antennas to be flagged.
        % Returns:
        %   img : structure containing generated images.
        %   sol : Solution structure as generated by the calibration routine.
        function [img, solx, soly] = genCalMap (obj, arrayconfig, flagant)
            load ('poslocal_outer.mat', 'poslocal');
            uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
            vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
			gparm.type = 'pillbox';
		    gparm.lim  = 0;
		    gparm.duv = 0.5;				% Default, reassigned from freq. of obs. to
											% image just the full Fov (-1<l<1)
		    gparm.Nuv = 512;				% size of gridded visibility matrix
		    gparm.uvpad = 512;				% specifies if any padding needs to be added
		    gparm.fft  = 1;

            img.tobs = obj.trecstart;
            img.freq = obj.freq;


	   	    solx = pelican_sunAteamsub (obj.acm_xx, img.tobs, img.freq, eye(obj.nelem), ... 
				flagant, obj.deb, 1, [], [], 'poslocal_outer.mat', [], []);
	   	    soly = pelican_sunAteamsub (obj.acm_yy, img.tobs, img.freq, eye(obj.nelem), ... 
				flagant, obj.deb, 1, [], [], 'poslocal_outer.mat', [], []);

            [uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
   		    [~, img.map_xx, ~, img.l, img.m] = ... 
		        fft_imager_sjw_radec (solx.calvis(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], img.tobs, img.freq, 0);
   		    [~, img.map_yy, ~, img.l, img.m] = ... 
		        fft_imager_sjw_radec (soly.calvis(:), uloc_flag(:), vloc_flag(:), ... 
					gparm, [], [], img.tobs, img.freq, 0);
        end;


        % Function to generate an uncalibrated map from the observed
        % visibilities. Returns generated images in the img structure
        % Both XX and YY images are returned.
        %
        % Arguments:
        %   arrayconfig : One of 'lba_outer', 'lba_inner'
        % Returns:
        %   img : structure containing generated images.
        function img = genUncalMap (obj, arrayconfig)
            assert (obj.freq ~= 0);
            load ('poslocal_outer.mat', 'poslocal');
            uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
            vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
			gparm.type = 'pillbox';
		    gparm.lim  = 0;
		    gparm.duv = 0.5;				% Default, reassigned from freq. of obs. to
											% image just the full Fov (-1<l<1)
		    gparm.Nuv = 512;				% size of gridded visibility matrix
		    gparm.uvpad = 512;				% specifies if any padding needs to be added
		    gparm.fft  = 1;

            img.tobs = obj.trecstart;
            img.freq = obj.freq;
   		    [~, img.map_xx, ~, img.l, img.m] = ... 
		        fft_imager_sjw_radec (obj.acm_xx(:), uloc(:), vloc(:), ... 
					gparm, [], [], img.tobs, img.freq, 0);
   		    [~, img.map_yy, ~, img.l, img.m] = ... 
		        fft_imager_sjw_radec (obj.acm_yy(:), uloc(:), vloc(:), ... 
					gparm, [], [], img.tobs, img.freq, 0);
        end;
			

	end % End of methods

end % End of classdef
