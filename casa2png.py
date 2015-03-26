# Script to convert casacore images generated by the AARTFAAC pipeline to FITS.
# pep/10Jun12

import glob;
import os;

print 'Generating PNG files from all CASAcore images in current folder';
default (imview)
 
os.system ('mkdir png');
# for image_file in glob.glob(os.path.join('fits/','*.fits')):
for imagename in glob.glob('*.image'):
	pngimage = imagename + '.png';
	print 'Converting file ', imagename, ' to ', pngimage;
	imview(raster=imagename, out=pngimage);
	print 'Done';
