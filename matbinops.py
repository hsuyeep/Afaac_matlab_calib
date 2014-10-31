""" Functions to write and read binary data in formats as generated by 
	various matlab scripts.
    pep/31Oct14
"""

import numpy as np;
import os;
import sys;
import struct;
import sys;
import time;
try:
	import pyfits;
	pyfitsFound = 1;
except ImportError:
	print '### pyfits not found!';
	pyfitsFound = 0;


class matBinImg:
	""" Abstraction of matlab written binary image """

	""" Read data written by wrimg2bin.m """
	def readimg2bin (self, fid):
		self._tobs  = np.fromfile(fid, dtype='float64', count=1);
		self._fobs  = np.fromfile(fid, dtype='float64', count=1);
		self._nlpix = np.fromfile(fid, dtype='float32', count=1);
		self._laxis = np.fromfile(fid, dtype='float32', count=self._nlpix);
		self._nmpix = np.fromfile(fid, dtype='float32', count=1);
		self._maxis = np.fromfile(fid, dtype='float32', count=self._nmpix);
		self._img   = np.fromfile(fid, dtype='float32', count=self._nmpix*self._nlpix);
		self._pol   = 'XX'; # NOTE: Hardcoded for now!
		self._ttype = 'mjdsec'; # Time is MJD, in seconds

	def getImg (self):
		self._img.shape = (self._nlpix, self._nmpix);
		return self._tobs, self._fobs, self._laxis, self._maxis, self._img;

#	default_fits_format_codes = {
#    np.bool_:'L', np.uint8:'B', np.int16:'I', np.int32:'J', np.int64:'K',
#    np.float32:'E', np.float64:'D', np.complex64:'C', np.complex128:'M'
#	}

	""" Returns the RA of a zenith pointing, given the observation time """
	def tobs2ra (self):
		if self._ttype.lower() == 'mjdsec':
			self._JD = self._tobs/86400. + 2400000.5; # MJDSec to JD
		elif self._ttype.lower() == 'jd':
            self._JD = self._tobs;
		elif self._ttype.lower() = 'utcdatenum':
            # Convert UTC time provided as datenum, to JD
            self._JD = self._tobs - datenum(1998, 2, 13, 12,0,0) + 2450858;
		else:
			print 'Unknown time type!';

		# Obtain LST for this time
		a1 = 24110.54841;
		a2 = 8640184.812866;
		a3 = 0.093104;
		a4 = -6.2e-6;
		polcoeff = [a4, a3, a2, a1];
		# Westerbork coordinates (From http://www.sonel.org/spip.php?page=gps&idStation=884)
		wsrt_lat = 52.91460037;
		wsrt_lon =  6.60449982;
		
		# Time in Julian centuries
		TU = (floor(JD) + 0.5 - 2451545) / 36525;
		
		# Greenwich Star Time in seconds
		GST = (JD - floor(JD) - 0.5) * 86400 * 1.002737811906 + polyval(polcoeff, TU);
		
		# Local Siderderial Time in radians
		# 240: number of seconds per degree
		LST =  ((GST + wsrt_lon*240) / 240) * pi / 180;
		# Convert to one cycle of 2*pi, which is 24 hrs;
		LST = LST - 2*pi*(floor (LST/(2*pi)));
		ra = LST - ha;



	
	def toFits(self, filename, clobber=False, axes=('ra--sin','dec--sin'),
		object='Zenith, Transit', telescope='AARTFAAC', 
		instrument='AARTFAAC', observer='', origin='AIPY',
		obs_date=time.strftime('%D'), cur_date=time.strftime('%D'), 
		ra=0, dec=0, d_ra=0, d_dec=0, epoch=2000.,
		freq=0, d_freq=0, bscale=0, bzero=0,history=''):
		"""Write image self._img to a FITS file.  Follows convention of VLA image
		headers.  "axes" describes dimensions of "self._img" provided.  (ra,dec) are
		the degree coordinates of image center in the specified "epoch". 
		(d_ra,d_dec) are approximate pixel-deltas for ra,dec (approximate because 
		if sine projections of these coordinates are specified--e.g. 
		"ra---sin"--then the deltas change away from the image center).  If a 
		"freq" axis is specified, then "freq" is the frequency of the first entry 
		(in Hz), and "d_freq" is the width of the channel.  The rest are pretty 
		self-explanitory/can be used however you want."""
	
		if pyfitsFound == 0:
			print '### Pyfits not found, unable to proceed.';
			return;
		
	
		self._img = self._img.squeeze()
		self._img.shape = (1,) * (len(axes) - len(self._img.shape)) + self._img.shape
		if len(self._img.shape) != len(axes):
			raise TypeError('to_fits: axes dimension list does not match self._img shape')
	
		phdu = pyfits.PrimaryHDU(self._img)
		phdu.update_header()
		phdu.header.update('TRANSPOS',0,comment='Import code for old AIPY convention.')
		phdu.header.update('OBJECT', object, comment='SOURCE NAME')
		phdu.header.update('TELESCOP', telescope)
		phdu.header.update('INSTRUME', instrument)
		phdu.header.update('OBSERVER', observer)
		phdu.header.update('DATE-OBS', obs_date,
		    comment='OBSERVATION START DATE DD/MM/YY')
		phdu.header.update('BSCALE ', bscale,
		    comment='REAL = FITS_VALUE * BSCALE + BZERO')
		phdu.header.update('BZERO  ', bzero)
		phdu.header.update('BUNIT  ', 'JY/BEAM ', comment='UNITS OF FLUX')
		phdu.header.update('EQUINOX', epoch, comment='EQUINOX OF RA DEC')
		phdu.header.update('EPOCH',epoch,'Epoch of coordinate system')
		phdu.header.update('DATAMAX', np.nanmax(self._img), comment='MAX PIXEL VALUE')# NOTE: Ignore nans
		phdu.header.update('DATAMIN', np.nanmin(self._img), comment='MIN PIXEL VALUE')# NOTE: Ignore nans
		print '--> Header population complete.';
	
		for i,ax in enumerate(axes):
			if ax.lower().startswith('ra') or ax.lower().startswith('glon'): val,delta = (ra, d_ra)
			elif ax.lower().startswith('dec') or ax.lower().startswith('glat'): val,delta = (dec, d_dec)
			elif ax.lower().startswith('freq'): val,delta = (freq, d_freq)
			elif ax.lower().startswith('stokes'): val,delta = (1, 1)
			else: val,delta = (0,0)
			phdu.header.update('CTYPE%d' % (i+1), ax.upper())
			if ax.lower().startswith('ra') or ax.lower().startswith('dec')\
			or ax.lower().startswith('glon') or ax.lower().startswith('glat'):
			    phdu.header.update('CRPIX%d' % (i+1), np.ceil(phdu.data.shape[-(i+1)]/2.))
			else:
				phdu.header.update('CRPIX%d' % (i+1), np.ceil(phdu.data.shape[-(i+1)]/2.))
			phdu.header.update('CRVAL%d' % (i+1), val)
			if (ax.lower().startswith('ra') or ax.lower().startswith('glon')) and delta>0:
				phdu.header.update('CDELT%d' % (i+1), delta*-1)
			else:
				phdu.header.update('CDELT%d' % (i+1), delta)
			phdu.header.update('CROTA%d' % (i+1), 0)
			phdu.header.update('NAXIS%d' % (i+1),phdu.data.shape[i])
	
		if history!='':
			history = [h.strip() for h in history.split("\n")]
			for line in history:
				if len(line)>1:
					if line.startswith('#'):
						for subline in word_wrap(line,80,1,1,'').split("\n"):
							phdu.header.add_history(subline)
					else:
						for subline in word_wrap(line,70,5,10,'#').split("\n"):
							phdu.header.add_history(subline)
		phdu.header.update('ORIGIN', origin);
		phdu.header.update('DATE', cur_date, comment='FILE WRITTEN ON DD/MM/YY');
		try:
			print '--> Writing to file ', filename;
			pyfits.writeto(filename,phdu.data,phdu.header,clobber=True, output_verify='warn');
		except:
			print '### Unable to write fits file to disk!';



class matBinVis:
	""" Abstraction of matlab written binary Covariance matrix """

	""" Read data written by wracm2bin.m, equivalent of readms2float.m """
	def readms2float (fid):
		return False;

# Test jig
if __name__ == '__main__':

	try: 
		import matplotlib.pylab as plt;
	except ImportError:
		print '### Matplotlib not found!';

	fid = open (sys.argv[1], 'rb');
	binimg = matBinImg ();

	binimg.readimg2bin (fid);
	tobs, fobs, l, m, img = binimg.getImg();
	print 'Time: %0.f, freq: %.0f' % (tobs, fobs);
	zendec = 52.5; # Latitude of AARTFAAC, good for zenith pointing
	zenra  = 23;  # Currently random, should come from time;
	ddec = -0.286480091781;
	dra  = ddec;

	# Try to write this out as a fits file.
	print 'Creating fits file...'
	try:
		binimg.toFits('testimage.fits', ra=zenra, dec=zendec, d_ra=dra, d_dec=ddec);
	except:
		print '### Unable to create fits file!';
	print 'Done';
	implt = plt.imshow (img);
	plt.colorbar();

#	while True:
#		binimg.readimg2bin (fid);
#		tobs, fobs, l, m, img = binimg.getImg();
#		print 'Time: %0.f, freq: %.0f' % (tobs, fobs);
#		implt.set_data(img);
#		plt.draw();
#		plt.pause (0.1);
