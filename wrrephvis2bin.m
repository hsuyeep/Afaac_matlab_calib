% Program to read a calibrated .bin file, rephase visibilities to the center 
% of a window of timeslices, and write out the rephased and calibrated 
% visibilities to file.
% Arguments:
% fname : Input .bin file containing uncalibrated visibilities
% offset: Offset in timeslices from the beginning of input file to start 
%		  calibrating.
%ntslices: Number of timeslices to calibrate
% NOTE: If output file is found to exist, the rephased visibilites are 
% appended to it.
% pep/03Dec12

function wrrephvis2bin (fname, offset, ntslices, tint)
    % Necessary definitions
    uvflag = eye (288); % Flag only the autocorrelations
	debuglev = 2;
	rephind = int32(tint/2);

	if ntslices == -1
		[ntslices, tmin, tmax, dt] = getnrecs (fname);
		disp (sprintf ('Found %d records at %f sec. resolution.', ntslices, dt));
	end
	% Generate output filename
	k = strfind (fname, '.bin');
	outfname = [fname(1:k-1) '_reph.bin'];
	disp (['Writing rephased visibilties to file :' outfname]);
		
	fin = fopen (fname, 'rb');
	fout = fopen (outfname, 'ab');
	[acct0, t, freq] = readms2float (fin, -1, -1);
	if (isempty(acct0) == true)
		disp ('wrrephvis2bin: Eof reached!');
		return;
	end;

	% Generate tint size containers.
	acc = zeros ([size(acct0), tint]);
    tobs = zeros (1, tint);    
	load ('poslocal.mat', 'poslocal','posITRF');

	for rec = 1:tint:ntslices
		disp (' ');
		disp (['###### Working on Rephased window :' num2str(rec)]);
		for ind = 1:tint
			[acct0, t, freq] = readms2float (fin, -1, -1);
			if (isempty(acct0) == true)
				disp ('wrrephvis2bin: Eof reached!');
				return;
			end;
			acc(:,:,ind) = acct0;
			tobs(ind) = t;
			disp (['Time: ' num2str(tobs(ind))]);
		end;

		for ind = 1:tint
			disp (['Rephasing from ' num2str(tobs(ind)) ' to ' ... 
					num2str(tobs(rephind))]);
			reph_acc = rephasetime(acc(:,:,ind), tobs(ind), tobs(rephind),... 
								freq, posITRF);
			reph_acc = single (reph_acc);


		% NOTE!NOTE!NOTE! THIS IS REQUIRED FOR wracm2bin TO FUNCTION CORRECTLY!
		% acc (antmask == 1) = NaN; % Flagged visibilities get a NaN.
		% whos acc;
			wracm2bin (fout, offset, reph_acc, [], tobs(ind), freq, -1);
		end
	end

	fclose (fin); fclose (fout);
 	if wrcalsol == 1
		fclose (fsol);
	end;
