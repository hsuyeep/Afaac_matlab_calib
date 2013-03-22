% Script to return a chosen timeslice from a .bin file generated by ms2weightfl.py
% Either an offset, or an absolute time can be given, with the function 
% returning the ACM, MJD time (in secs), Freq (in Hz), and the weight 
% associated with the returned ACM. This weight is generated by factoring in
% data losses to the correlator, as well as the visibilities lost due to 
% flagging over an integration time period.
% NOTE: If recoffset=-1, implies return the next record, for any other 
% positive number, return the record with that number offset from beginning
% of file.

% Arguments:
%   fid      : fid of the binary data file.
%   recoffset: Offset (in units of records) of the desired timeslice.
%   rettime  : Time in MJD seconds of the desired timeslice.
% Returns  :
%   acc      : The complex, Hermitean symmetric Array Correlation Matrix with 
%              autocorrelations removed,corresponding to chosen timeslice.
%  weight	 : Weight associated with every visibility in the acc, as a square
%				matrix, including the autocorrelations.
%   tobs     : Time of observation, in MJD seconds, as a double.
%   freq     : Frequency of observation, in Hz.
% pep/01May12
% function [acc, tobs, freqobs] = readms2float (filename, recoffset, rettime)

function [acc, weight, tobs, freqobs] = readms2weightfl(fid, recoffset, rettime, Nelem)
	if (isempty(Nelem) == 1) 
		Nelem = 288; 
	end;

	nblines = Nelem * (Nelem + 1)/2; 
	tobs = 1; freqobs = 1;
	recsize = 8*(2+nblines); % Bytes, assumed double =8, float=4 bytes.
	% fid = fopen (filename, 'rb');

	if (recoffset > 0)
		fseek (fid, (recoffset-1)*recsize, 'bof');
		tobs = fread (fid, 1, 'double');
		if (isempty (tobs) == true)
 			disp ('readms2float: EoF reached!'); 
			acc = []; tobs = []; freqobs = [];
			return;
		end;
	    freqobs = fread (fid, 1, 'double');
        % Reading real and imaginary, available as a stream of floats.
        % even floats being real parts, odd floats being imag
		a  = fread (fid, 2*nblines, 'float'); % Read one ccm worth
		disp (['Time at offset ', num2str(recoffset), ' recs: ', ...
			  num2str(tobs)]);
	elseif (recoffset == -1) % Return next record
		tobs = fread (fid, 1, 'double');
		if (isempty (tobs) == true)
 			disp ('readms2float: EoF reached!'); 
			acc = []; tobs = []; freqobs = [];
			return;
		end;
	    freqobs = fread (fid, 1, 'double');
		a  = fread (fid, 2*nblines, 'float'); % Read one ccm worth.
		w = fread (fid, nblines, 'float'); % Read in the weights.
		% disp (['Time: ', num2str(tobs), ' Freq: ', num2str(freqobs)]);
	end	

    comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex

    % create instantaneous ccm from vector
    acm = triu (ones (Nelem));
	weight = triu (ones (Nelem));
    acm (acm == 1) = comp;
	weight(weight == 1) = w;

    % Removing autocorrelations
    acc = acm + acm' - diag(diag(acm));
	%fclose (fid);
end