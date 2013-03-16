% Script to return the number of records in a binary ACM file.
% pep/18Oct12

function [ntimes, tmin, tmax, dt] = getnrecs (fname)
	D = dir (fname);
    fid = fopen (fname, 'rb');
	fseek (fid, 0, 'bof'); % Go to beginning of file
	Nelem = 288; 
    nblines = Nelem * (Nelem + 1)/2;
    tobs = 1; freqobs = 1;

    recsize = 8*(2+nblines); % NOTE: Works only for single channel data!
	ntimes = D.bytes / recsize;
	[acc, tmin, freq] = readms2float (fid, -1, -1, 288);
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	fseek (fid, -recsize, 'eof'); % Go to the end of the file.
	[acc, tmax, freq] = readms2float (fid, -1, -1, 288);
	dt = tobs - tmin;
	fclose (fid);


