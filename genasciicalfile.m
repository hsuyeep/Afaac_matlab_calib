% Script to generate a calibration solution file from a calsol.bin file
% Used in the one time caibration of the QuickLook images.
% pep/26Jun15
% Arguments:
%  calfile: file name of the calibration solution binary file.
%  toff   : record offset to read.
% Returns:
%  rec : calsol record.

function [rec] = genasciicalfile (fname, pol, array, toff)
	fid = fopen (fname, 'rb');
    if fid < 0
		error ('genasciicalfile: Error in opening file.\n');
	end;

	rec = readcalsol (fid);	
	tstamp = datestr (mjdsec2datenum(rec.tobs), 30);
	calfname = [tstamp, '_', pol, '_', array, '.cal'];
	fout = fopen (calfname, 'wt');
	fprintf (fout, '# Time: %s, Freq: %.2f, Pol: %s, Array: %s\n', tstamp, rec.freq, pol, array);
	fprintf (fout, '# Antenna  Gain_real  Gain_imag\n');
	for ind = 1:length (rec.gainsol)
		fprintf (fout, '%03d  %7.2f %7.2f\n', ind, real(rec.gainsol(ind)), imag(rec.gainsol(ind))); 
	end;
