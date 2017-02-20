""" Script to write images generated in local topocentric coordinates
    as FITS files, with a changing pointing center (zenith) based on 
	observation time.
"""

import sys;
import numpy as np;
import aipy as a;
import pyfits;
import matplotlib.pylab as plt;
import matbinops;

if __name__ == '__main__':
	fid = open (sys.argv[1], 'rb');
	skip = int (sys.argv[2]);

	print '<-- Skipping ', sys.argv[2], ' secs.';
	# We assume image files have been written by writeimg2bin.m
	binimg = matbinops.matBinImg ();

	while True:
		for ind in range (0, 10):
			binimg.readimg2bin (fid);
		tobs, fobs, l, m, img = binimg.getImg();
		print 'Time: %0.f, freq: %.0f' % (tobs, fobs);
		binimg.toFits ();
