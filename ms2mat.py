#!/usr/bin/python
# Program to extract a subset of visibility data from .MS files, 
# and write them out to multiple .mat files, one per timeslice.
# Meta data (antenna names, spectral windows, bandwidths etc. are 
# written out into another .mat file, and is common for the full dataset.
# pep/03Jul15

import traceback
try:
	import sys;
	import argparse;
	from pyrap.tables import *
	from numpy import *
	from scipy.io import *
except ImportError:
	print 'Pyrap module not found!';

import numpy;
import os;

def main ():
#	if (len(sys.argv) < 2):
#		print 'Convert a .MS into a series of .mat files. Meta information available separately.';
#		print 'Usage: ', sys.argv[0], ' filename.MS [rectimesize]';
#		print '	rectimesize = number of timeslices in a single matlab variable (default = 1), UNIMPLEMENTED!';
#		sys.exit (-1);

	parser = argparse.ArgumentParser()
	parser.add_argument("--ms", help="Specify the MeasurementSet to be converted.");
	parser.add_argument("--t2file", help="Specify the number of timeslices to write to a single .mat file (Currently UNIMPLEMENTED). Default is a single timeslice", default='1');
	parser.add_argument("--trng", help="Specify the time range (as a python list) to extract from the MS, e.g., as MJDstart:MJDend (Currently UNIMPLEMENTED). Default is all timeslices.", default='-1');
	parser.add_argument("--frng", help="Specify the spw range (as a python list) to extract from the MS, e.g., as as spwstart:spwend (Currently UNIMPLEMENTED). Default is all spw.", default='-1');
	parser.add_argument("--instr", help="Specify the source instrument of this MS. Default is all 'APERTIF'.", default='APERTIF');
	client = parser.parse_args();

	# Get rid of trailing '/' to allow correct basename.
	if (client.ms[-1] == '/'):
		client.ms = client.ms[:-1];
	print '\n<-- Operating on file ', client.ms;
	tab = table (client.ms); 

	tblksize = int(client.t2file);
	
	# Get total number of antennas in the antenna table.
	anttab = table (tab.getkeyword ('ANTENNA'));
	nelem  = len (anttab.getcol ('NAME'));
	specttab = table (tab.getkeyword ('SPECTRAL_WINDOW'));
	nspw   = len(specttab.getcol ('NAME'));
	nbline = nelem * (nelem+1)/2;  # Including autocorrelations.
	nrec   = len (tab.getcol ('TIME')); # = ntimes*npol*nchan;
	dt     = tab.getcol ('EXPOSURE')[0];
	ntslot = int (round ((tab [nrec-1]['TIME'] - tab [0]['TIME'])/dt)) + 1; 

	if (nrec != (ntslot*nbline*nspw)):
		print '### Total records not consistent with total timeslices and total baselines!';
		print 'ntslots: %d, nrec :%d, ntslot*nbline=%d'% (ntslot, nrec, ntslot*nbline);

	npol   = 4; # Currently hardcoded!
	nchan  = specttab.getcol ('NUM_CHAN')[0];
	print '    Found ', nelem, ' antennas in MS.';
	print '    Found ', nspw, ' Spectral Windows in MS.';
	print '    Time range:',tab [0]['TIME'], '-', tab [nrec-1]['TIME'], '(', ntslot, ') secs';
	print 'Exposure  :', tab.getcol('EXPOSURE')[0],' secs';

	for ind, spw in enumerate(specttab.getcol ('NAME')):
		print '<-- SPW: %s, Total band: %.2f Hz, central freq: %.2f Hz, Channels: %d.'% (spw, specttab.getcol('TOTAL_BANDWIDTH')[ind],specttab.getcol('REF_FREQUENCY')[ind], nchan);
	print '<-- Splitting ', ntslot, 'timeslices into blocks of ', tblksize;

	matfilename = os.path.basename (client.ms).split('.')[0] + '_meta.mat';
	print '<-- Writing out one-time meta data to matlab file: ', matfilename;
	# ITRF postions of reference within station. NOTE: Antenna phase ref. is same as station position.
	antpos = zeros ( (nelem, 3), dtype=double); 
	antname = [];
	for ant in range (0, nelem):
		antname.append (anttab [ant]['NAME']);
		antpos [ant] = anttab [ant]['POSITION'];

	chan2freq = zeros ( (nspw, nchan), dtype=double);
	chanbw = zeros ( (nspw, nchan), dtype=double);
	for ind,spw in enumerate (specttab.getcol('NAME')):
		chan2freq[ind] = specttab [ind]['CHAN_FREQ'];
		chanbw[ind] = specttab [ind]['CHAN_WIDTH'];

	# Create data structures for each timeslice
	acm = zeros ( (tblksize, nspw, nelem, nelem, nchan, npol), dtype=complex64); # The ACM
	uvw = zeros ( (tblksize, nspw, nelem, nelem, 3), dtype=double);        # UVW coordinates of each visibility
	lofar_fullres_fl = zeros ( (1, 8), dtype=uint16);      # LOFAR_FULL_RES_FLAG

	rec = dict ();
	rec ['ANTPOS'] = antpos;
	rec ['ANTNAME'] = antname;
	rec ['SPW_CHAN_FREQ'] = chan2freq;
	rec ['SPW_CHAN_BW'] = chanbw;
	rec ['SPW_NAME'] = specttab.getcol('NAME');
	rec ['DATA_DIM'] = ['NTIMESLOTS', 'NSPW', 'NANT', 'NANT', 'NCHAN', 'NPOL'];
	rec ['UVW_DIM'] = ['NTIMESLOTS', 'NSPW', 'NANT', 'NANT', 'U', 'V', 'W'];
	rec ['SPW_CHAN_FREQ_DIM'] = ['NSPW', 'NCHAN'];
	rec ['SPW_CHAN_BW_DIM'] = ['NSPW', 'NCHAN'];
	savemat (matfilename, rec);
	
	for tblk in range (0, ntslot/tblksize):

		# Create a single matlab variable per timeslice
		tacm = tab[tblk*nbline*nspw]['TIME'];
		matfilename = os.path.basename (client.ms).split('.')[0] + '_' + str(int(tacm))+'.mat';
		print '<-- Now writing timeslot %d to file %s.' % (tacm, matfilename); 
		for spw in range (0, nspw):
			for vis in range (0, nbline):
				ind = tblk * nbline * spw + vis;
				a0 = tab [ind]['ANTENNA1']; 
				a1 = tab [ind]['ANTENNA2'];	
				trec = tab [ind]['TIME'];
		
				acm [0, spw, a0, a1] = tab [ind]['DATA'];
				acm [0, spw, a1, a0] = acm [0, spw, a0, a1].conjugate();
				uvw [0, spw, a0, a1] = tab [ind]['UVW'];
				uvw [0, spw, a1, a0] = -uvw [0, spw, a0, a1];
				if (client.instr == 'LOFAR'):	
					lofar_fullres_fl = tab [ind]['LOFAR_FULL_RES_FLAG'];
		
		rec = dict();
		rec ['TIME'] = tacm;
		rec ['DATA'] = acm;
		rec ['UVW'] = uvw;
		if (client.instr == 'LOFAR'):	
			rec ['FULLRESFL'] = lofar_fullres_fl;
	
		savemat (matfilename, rec);
	
	return;

if __name__ == "__main__":
	main ();
