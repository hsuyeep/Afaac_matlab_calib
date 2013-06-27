#!/usr/bin/python
# Program to extract complete visibility data from .MS files, 
# and write them out to multiple .mat files, one per timeslice.
#  The frequency and time information is written out to a
# separate .mat file.
# Output format:
#   Each .mat file contains the following variables:
#	'DATAXX': Nant x Nant x Nchan
#	'DATAXY': Nant x Nant x Nchan
#	'DATAYX': Nant x Nant x Nchan
#	'DATAYY': Nant x Nant x Nchan
#	'UVW'   : Nant x Nant x 3 

# pep/26Jun13

import traceback
try:
	import sys;
	from pyrap.tables import *
	from numpy import *
	from scipy.io import *
except ImportError:
	print 'Pyrap module not found!';

import numpy;
import os;

def main ():
	if (len(sys.argv) < 2):
		print 'Usage: ', sys.argv[0], ' filename.MS [rectimesize]';
		print '	rectimesize = number of timeslices in a single matlab variable (default = 1), UNIMPLEMENTED!';
		sys.exit (-1);

	msfile = sys.argv[1]; 
	tab = table (sys.argv[1]); 

	tblksize = 1; # Default
	if (len(sys.argv) == 3):
		tblksize = int(sys.argv[2]);
	

	# Get total number of antennas in the antenna table.
	anttab = table (tab.getkeyword ('ANTENNA'));
	nelem  = len (anttab.getcol ('NAME'));
	specttab = table (tab.getkeyword ('SPECTRAL_WINDOW'));
	nbline = nelem * (nelem+1)/2; 
	nrec   = len (tab.getcol ('TIME')); # = ntimes*npol*nchan;
	dt     = tab.getcol ('EXPOSURE')[0];
	ntslot = int ((tab [nrec-1]['TIME'] - tab [0]['TIME'])/dt); 

	if (nrec != (ntslot*nbline)):
		print '### Total records not consistent with total timeslices and total baselines!';
		print 'ntslots, nrec ;', ntslot, ',', nrec, 'ntslot*nbline=', ntslot*nbline;

	npol   = 4; # Currently hardcoded!
	nchan  = specttab.getcol ('NUM_CHAN')[0];
	print '\n Found ', nelem, ' antennas in MS.';
	print '\n\nTime range:',tab [0]['TIME'], '-', tab [nrec-1]['TIME'], '(', ntslot, ') secs';
	print 'Exposure  :', tab.getcol('EXPOSURE')[0],' secs';
	print 'Total band:', specttab.getcol('TOTAL_BANDWIDTH')[0],' Hz, around ',specttab.getcol('REF_FREQUENCY')[0], ' Hz.';
	print 'Subband   :', specttab[0]['NAME'];
	print 'Channels  :', nchan;
	print 'Splitting ', ntslot, 'timeslices into blocks of ', tblksize;

	print 'Writing out one-time meta data:'
	# ITRF postions of reference within station. NOTE: Antenna phase ref. is same as station position.
	antpos = zeros ( (nelem, 3), dtype=double); 
	antname = [];
	for ant in range (0, nelem):
		antname.append (anttab [ant]['NAME']);
		antpos [ant] = anttab [ant]['POSITION'];

	chan2freq = zeros ( (1, nchan), dtype=double);
	chanbw = zeros ( (1, nchan), dtype=double);
	chan2freq = specttab [0]['CHAN_FREQ'];
	chanbw = specttab [0]['CHAN_WIDTH'];
	rec = dict ();
	rec ['ANTPOS'] = antpos;
	rec ['ANTNAME'] = antname;
	rec ['CHAN_FREQ'] = chan2freq;
	rec ['CHAN_BW'] = chanbw;
	savemat (specttab[0]['NAME'] + '_meta.mat', rec);
	
	# Create data structures for each timeslice
	acm = zeros ( (nelem, nelem, nchan, npol), dtype=complex64); # The ACM
	uvw = zeros ( (nelem, nelem, 3), dtype=double);        # UVW coordinates of each visibility
	lofar_fullres_fl = zeros ( (1, 8), dtype=uint16);      # LOFAR_FULL_RES_FLAG

	# Get all relevant information: Can we read in just a timeslice worth?
	# trow = tab.row (['TIME', 'DATA', 'UVW', 'ANTENNA1', 'ANTENNA2']); 

	for tblk in range (0, ntslot/tblksize):
	# for tblk in range (0,5): 

		# Create a single matlab variable per timeslice
		tacm = tab[tblk*nbline]['TIME'];
		print 'Timeslot ', tacm;
		for vis in range (0, nbline):
			ind = tblk * nbline + vis;
			a0 = tab [ind]['ANTENNA1']; 
			a1 = tab [ind]['ANTENNA2'];	
			trec = tab [ind]['TIME'];
	
			if trec != tacm:
				print 'Misalignment! Encountered a new time within ACM! tacm = ', tacm, 'trec = ', trec;
	
			acm [a0, a1] = tab [ind]['DATA'];
			acm [a1, a0] = acm [a0, a1].conjugate();
			uvw [a0, a1] = tab [ind]['UVW'];
			uvw [a1, a0] = -uvw [a0, a1];
			lofar_fullres_fl = tab [ind]['LOFAR_FULL_RES_FLAG'];
		
		rec = dict();
		rec ['TIME'] = tacm;
		rec ['DATA'] = acm;
		rec ['UVW'] = uvw;
		rec ['FULLRESFL'] = lofar_fullres_fl;
	
		savemat (specttab[0]['NAME'] + '_' + str(int(tacm)), rec);
	
	return;

if __name__ == "__main__":
	main ();
