% Program to apply an average cal. solution to uncalibrated visibilities, based
% on instantaneous gain solutions. Generates a calibrated .bin file.
%
% Arguments:
% 	fname : Input .bin file containing uncalibrated visibilities
% calsolfname: Input .bin file containing calibration solutions.
% avgover : Timeslices over which to average calib. solutions, before applying.
% 	offset: Offset in timeslices from the beginning of input file to start 
%		  	calibrating. (NOTE: Calsol file should also have the same time 
%			resolution as the uncalibrated visibilities!)
% ntslices: Number of timeslices on which to apply calibration solutions.
% peelon  : Bool controlling the turning on of Peeling. (TODO)
% appavg  : Bool controlling the application of averaged solutions to
% 			uncalibrated visibilities and their writing out to file. If 0, only
% 			averaging is carried out.
% usrflagant : Antennas to flag, additional to the ones in the calibration
%           solution per timeslice.
% Returns:
%	 avg_re_gain/avg_im_gain: The mean gain, averaged over the averaging window.
% pep/20Jan13

function [gain_acc] = applyavgcal2bin (fname, calsolfname, ... 
								avgover, offset, ntslices, peelon, appavg, usrflagant)


    % Necessary definitions
	Nelem = 288;
    uvflag = eye (Nelem); % Flag only the autocorrelations
	nblines = Nelem * (Nelem + 1)/2; 

    % flagant = [6, 103];  % For day-time observations over 3hr (2011 data)

	% For LBA_INNER_BAND60 data
    % flagant = [1:12, 47, 48, 95,96, 143, 144, 191, 192, 239, 240, 287, 288 ];   
    % flagant = [1:12, 51, 193, 239,273 ]; % For LBA_OUTER_BAND60 data
    % flagant = [1:48, 51, 239]; % For LBA_OUTER_BAND_SPREAD data
    % flagant = [49:96, 239, 241:288]; % for 03285_dawn_spread data
	% disp ('NOTE! NOTE: Currently, the list of antennas to be flagged is a parameter within the script, and needs to be changed for different datasets!');
	% disp (['Currently flagged antennas: ' flagant]);
	
	debuglev = 0;

	if ntslices == -1
		[ntslices, tmin, tmax, dt] = getnrecs (fname);
		fprintf (1, 'applyavgcal: Found %d uncalibrated records at %f sec. resolution.\n',...
				ntslices, dt);
	end

	% File operations.
	fin = fopen (fname, 'rb');
	if (fin < 1)
		fprintf (2, '### Error in opening %s.\n', fname);
		return;
	end;
	fsol = fopen (calsolfname, 'rb');
	if (fsol < 1)
		fprintf (2, '### Error in opening %s.\n', calsolfname);
		return;
	end;

	% Generate output filename, operate on various files.
	if (appavg == 1)
		k = strfind (fname, '.bin');
		outfname = sprintf ('%s_%02davgcal.bin', fname(1:k-1), avgover); 
		fprintf (1, '<-- applyavgcal: Writing Calibrated visibilties to file %s.\n', ...
				 outfname);
		fout = fopen (outfname, 'wb');
		if (fout < 1)
			fprintf (2, '### Error in opening file %s\n.', outfname);
			return;
        else
            fprintf (1, '<-- Writing calibrated visibilities to file %s\n', outfname);
		end;
	end;

	% Calculate visibility record size after reading in one record
	[acc, tobs, freq] = readms2float (fin, offset, -1, 288);
	visrecsize = 8*(2+nblines); % Bytes, assumed double =8, float=4 bytes.
	if (isempty(acc) == true)
		disp ('applyavgcal: Eof reached!');
		return;
	end;
	acc = single (acc);

	% Read in a gainsol record, generate datastructures
	rec = readcalsol (fsol);
	if (isempty (rec.tobs) == 1)
		disp('applyavgcal2bin: End of solution file reached!'); 
		return;
	end;
    
   	antmask = zeros (size (acc));
	gainmask = zeros (1, size (acc,1));
    % rem_ants = length(acc) - length(rec.flagant);
	recsize = whos ('rec');
	ncycles = int32(ntslices/avgover);
    fprintf (1, '<-- Avraging over %d timeslices, found %d cycles.\n', avgover, ncycles);
    gain_acc = zeros (Nelem, ncycles)+i*zeros(Nelem, ncycles);
	
	wr_acm = zeros (size (acc, 1), size (acc, 2));

	% Check calibration solution file for timeslices corresponding to visibility files.
	% Move to matching timeslice, if required.
	fprintf (1, '<-- Vis. Time: %f, Gainsol. Time: %f', tobs, rec.tobs);
	toffset = int32 (tobs - rec.tobs); % Round to the nearest record.
	if (toffset < 0) % Have to move in visibilities!
		if (fseek (fin, abs(toffset)*visrecsize, 'cof') < 0)
			 [msg, num] = ferror (fin); disp (['Seek error! ' msg]);
		end;
	elseif (toffset > 0) % Have to move in solutions!
		if (fseek (fsol, toffset*recsize.bytes, 'cof') < 0 )
			 [msg, num] = ferror (fsol); disp (['Seek error! ' msg]);
		end;
	end;
	
	% Check for matching number of antennas in visibilities and gainsolutions.
	% TODO: Currently we work assuming the same number of flagged antennas
	% in the visibility set and the gain solutions.
	% if (length(flagant) > rec.gainsol_len) % We have more solutions than required!
		% Set 
	% end;

	% Read in new cal. solution
	rec = readcalsol (fsol);

	% Read in new uncalibrated visibilities
	[acc, tobs, freq] = readms2float (fin, offset, -1, 288);
	
	disp (sprintf('After align: vis. time: %f, sol. time: %f', tobs, rec.tobs));

    flagant = union (usrflagant, rec.flagant);

	% Create a flagged mask for dealing with flagged antennas.
	disp (['## Flagging dipole numbers : ' num2str(flagant')]);
    for ind = 1:length(flagant)
    	antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
	end
	gainmask (flagant) = 1;
	rem_ants = length(acc) - length(flagant);

	% Operate on input visibilities.
    sigman_acc = zeros ([size(rec.sigman) ncycles]);

	for avgcycle = 1:ncycles
		% Average gainsolution parameters over desired timeslices.
		% TODO: Do sigma clipping over timeslices.
		for ind = 1:avgover
			try
				rec = readcalsol (fsol);
			catch err
				fprintf (2, 'EoF reached!');
				break;
			end;

            gain_acc (:,avgcycle) = gain_acc (:,avgcycle) + rec.gainsol;
            sigman_acc (:,:,avgcycle) = sigman_acc (:,:,avgcycle) + rec.sigman;
		end
        gain_acc (:,avgcycle) = gain_acc (:,avgcycle) ./ avgover;
        sigman_acc(:,:,avgcycle) = sigman_acc (:,:,avgcycle) ./ avgover;
	 % try
		% Apply averaged calibration solutions to uncalibrated visibilities
		if (appavg == 1)
			for ts = 1:avgover
				disp (' ');
				disp (['###### Working on Timeslice ' num2str((avgcycle-1)*avgover+ts+offset)]);
				[acc, tobs, freq] = readms2float (fin, -1, -1, 288);
				if (isempty(acc) == true)
					disp ('wrcalvis2bin: Eof reached!');
					return;
				end;
	    		acc = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
				acc = single (acc);
				disp (['Time: ' num2str(tobs)]);
	
				% Calibrate!
				% calvis = (gain_acc(gainmask ==0,avgcycle) * gain_acc(gainmask==0,avgcycle)') .* (acc -  sigman_acc(:,:,avgcycle));
                calvis = (gain_acc(gainmask ==0,avgcycle) * gain_acc(gainmask==0,avgcycle)') .* (acc);                
                wr_acm (antmask == 0) = single (calvis); 
				% calvis = single (calvis);
	
	
				%NOTE!NOTE!NOTE!THIS IS REQUIRED FOR wracm2bin TO FUNCTION CORRECTLY!
				% wr_acm (antmask == 1) = NaN; % Flagged visibilities get a NaN.
				wr_acm (antmask == 1) = 0; % Flagged visibilities get a 0, as NaN screws up svds...
				wracm2bin (fout, wr_acm, flagant, tobs, freq);
			end;
		end;

     % catch excep
%		disp ('Exception caught! Clean quitting...');
%		getReport (excep, 'basic');
%	    % fclose (fout); fclose (fin);	
%	 end;
	end;

	if (fin  > 0) fclose (fin);  end;
	if (appavg == 1 && fout > 0) fclose (fout);  end;
	if (fsol > 0) fclose (fsol);  end;
