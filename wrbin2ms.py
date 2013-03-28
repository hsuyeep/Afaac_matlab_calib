#!/usr/bin/python
# Program to accept a .bin and a .MS file as input, and to substitute the 
# visibilities in the MS file with the visibilities in the .bin file.
# NOTE: THE ORIGINAL MS FILE IS UPDATED! ORIGINAL IS LOST!
# This is useful for generating calibrated visibilites using SJW's Matlab code, 
# while the imaging is carried out in CASA.
# pep/18Oct12

import traceback
try:
	from pyrap.tables import *
except ImportError:
	print 'Pyrap module not found!';

import numpy;
import os;
import sys;
import struct;


def main ():
	if (len(sys.argv) < 3):
		print 'Usage: ', sys.argv[0], ' filename.bin filename.MS';
		sys.exit (-1);

	msfile = sys.argv[1]; 
	binfile = sys.argv[2]; # Binary file containing calibrated visibilities.

	nelem  = 288;        # TODO: Get from MS
	nblines= nelem*(nelem+1)/2; 
	nsubs  = 1;          # NOTE: Currently catering only to a single subband!
						 # Number of subbands
	
	# Bin file manipulation
	recsize = 8*(2+nbl); # rec size in bytes.
    siz = os.stat (binfile).st_size;
    nrec = siz/recsize;
    print 'Binary File: ', binfile, '\nSize: ', siz, 'Bytes. Records: ', nrec;
    fbin = open (binfile, 'rb');
	

	# MS manipulation.
	# Create an array for every timeslice from all subbands
	tobs = numpy.zeros (1, 'd'); # Time of obs, as double
	acm = numpy.zeros (nsubs*2*nblines, 'f');
	# TODO: Dont try reading the entire MS into mem! Go in chunks
	tsb0 = taql ('select DATA from $msfile');
	ttim0 = taql ('select TIME from $msfile');
	tab = table (sys.argv[1]); # Need to open table for getting spectral info
	tab1 = table(tab.getkeyword('SPECTRAL_WINDOW'));
	nchan = tab1[0]['NUM_CHAN'];
	print '-->Found ', nchan, ' Spectral channels';
	freqobs = numpy.zeros (nchan, 'd'); # Freq. of obs, Hz, as double
	for i in range (0, nchan):
		freqobs[i] = tab1[0]['CHAN_FREQ'][i];

	# Dont need the tables anymore
	# tab1.close(); tab.close()

	if nchan >  1:
		print '### Currently handle only 1 channel!';
		sys.exit(-1);

	ntimes = tsb0.nrows()/nblines;

	# Read in the first binary record.
	rec = fbin.read (recsize);
    tim = struct.unpack ('d', rec[0:8]);   # Returns a tuple
    freq = struct.unpack ('d', rec[8:16]); # Returns a tuple

	print '---> Ants = ',nelem,' Subbands = ',nsubs,' Times = ', ntimes;

	print 'Time of first record in MS : %10.1f' % (ttim0[0]['TIME']);
	print 'Time of first record in bin: %10.1f' % (tbin[0]);


#	i = 0;
#	done = 0;
#	for i in range (0, ntimes):
#		if done == 1:
#			print 'Done: ',i, 'timeslices' ;
#			break;
#	
#	
#		print 'Time:%f Freq:%f' % (ttim0[i*nblines+1]['TIME'],freqobs[0]) ;
#		tobs[0] = ttim0[i*nblines + 1]['TIME'];
#		for j in range (0, nblines): # NOTE: Previously had an off-by-one error! pep/27Apr12
#			try:
#				acm[2*j  ]=tsb0[i*nblines+j]['DATA'][0][0].real;
#				acm[2*j+1]=tsb0[i*nblines+j]['DATA'][0][0].imag;
#				# print 'j:', '%05d' % j, '(','%8.4f' % tsb0[i*nblines+j]['DATA'][0][0].real,',', '%8.4f' % tsb0[i*nblines+j]['DATA'][0][0].imag, ')';
#			except KeyboardInterrupt:
#				print '    ###  Keyboard interrupt...';
#				print '    ### Quit after current timeslice!';
#				print '(i,j)=',i,j;
#				done = 1;
#
#		tobs.tofile (ffloat);
#		freqobs.tofile(ffloat);
#		acm.tofile  (ffloat);
#		bytes = tobs.nbytes + freqobs.nbytes + acm.nbytes;
	
	return;

if __name__ == "__main__":
	main ();
