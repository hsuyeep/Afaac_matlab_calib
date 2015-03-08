#! /usr/bin/env python
"""
This is a general-purpose script for making faceted, spherical maps (stored in
Healpix FITS format) from individual "flat" maps stored in FITS files.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse, ephem

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.fits')
o.set_description(__doc__)
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created. If multiple frequnecies being written, append ID to end of filename')
o.add_option('-n', '--nobeam', dest='nobeam', action='store_true',
    help='Do not apply a Gaussian downweighting')
o.add_option('--min_alt', dest='min_alt', default=20.0, type='float',
    help='Mask out pixels below a minimum altitude (degrees), default=20.0')
o.add_option('--nside', dest='nside', type='int', default=256,
    help='NSIDE parameter for map, if creating a new file. default:256')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
o.add_option('-S','--stokes', dest='stokesid', default='I', type='str',
    help='Stokes Parameter to use (I,Q,U,V), default=I')
o.add_option('-F','--freq', dest='freqid', default='0', type='str',
    help='Frequency slice to use, all: produce individual maps from each frequency channel, default=0')
opts, args = o.parse_args(sys.argv[1:])

stokesDict={'i':0,'q':1,'u':2,'v':3}
stokesId=stokesDict[opts.stokesid.lower()]
freqId=opts.freqid

prev_dra = None
for i, filename in enumerate(args):
    img0, kwds = a.img.from_fits(filename)
    if freqId.lower().startswith('a'):
        freqIdList=range(img.shape[0])
    else:
        freqIdList=map(int,freqId.split(','))
    for fid in freqIdList:
        # Open skymap
        if len(freqIdList)>1: mapFn=opts.map.split('.hpx')[0]+'.f%i.hpx'%fid
        else: mapFn=opts.map
        if os.path.exists(mapFn): skymap = a.map.Map(fromfits=mapFn)
        else: skymap = a.map.Map(nside=opts.nside)
        skymap.set_interpol(opts.interpolate)

        # import pdb; pdb.set_trace();

        #img=img[freqId,stokesId]
        # img=img0[fid,stokesId] # Removed by Peeyush for image files without freq/stokes axis.
        img = img0;
        print 'Shape of img0: ' , n.shape (img0), ' img: ',  n.shape (img);
        img = img.squeeze()
        print '2. Shape of img: ', n.shape (img);
        #transpose and flip image to correct direction
        #img=n.fliplr(n.flipud(img.T))
        img=n.fliplr(n.flipud(img))
        print 'MIN',n.min(img),n.max(img)
        # Read ra/dec of image center, which are stored in J2000
        assert(kwds['epoch'] == 2000)
        s = ephem.Equatorial(kwds['ra']*a.img.deg2rad, kwds['dec']*a.img.deg2rad, 
            epoch=ephem.J2000)
        # To precess the entire image to J2000, we actually need to precess the
        # coords of the center back to the epoch of the observation (obs_date),
        # get the pixel coordinates in that epoch, and then precess the coords of
        # each pixel.  This is because extrapolating pixel coords from the J2000 
        # center assumes the J2000 ra/dec axes, which may be tilted relative to
        # the ra/dec axes of the epoch of the image.
        obs_date=' '.join(kwds['obs_date'].split('T'))
        ##HACK: for some reason the FITS header has the wrong time, so for now we get the
        ##JD from the file name and convert that to an obs time
        #ifn=filename.split('/')[-1]
        #ts=float(ifn.split('.')[1]+'.'+ifn.split('.')[2])
        #jd=ts-2415020
        #correct_obs_date=ephem.Date(jd)
        #obs_date=correct_obs_date

        s = ephem.Equatorial(s, epoch=obs_date)
        ra, dec = s.get()
        print '-----------------------------------------------------------'
        print 'Reading file %s (%d / %d)' % (filename, i + 1, len(args))
        print kwds
        print 'Pointing (ra, dec):', ra, dec
        print 'Image Power:', n.abs(img).sum()
        if prev_dra != kwds['d_ra']:
            prev_dra = kwds['d_ra']
            DIM = img.shape[0]
            RES = 1. / (kwds['d_ra'] * a.img.deg2rad * DIM)
            im = a.img.Img(DIM*RES, RES, mf_order=0)
            tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
            # Define a weighting for gridding data into the skymap
            map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
            if opts.nobeam: map_wgts=n.ones_like(map_wgts)
            #add a mask to pixels below minimum altitude
            min_alt=opts.min_alt    #degrees
            mask_ind=n.where(tz < n.sin(min_alt*a.img.deg2rad))
            map_wgts.mask[mask_ind]=True
            map_wgts.shape = (map_wgts.size,)
            valid = n.logical_not(map_wgts.mask)
            map_wgts = map_wgts.compress(valid)
        # Get coordinates of image pixels in the epoch of the observation
        ex,ey,ez = im.get_eq(ra, dec, center=(DIM/2,DIM/2))
        ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
        img = img.flatten(); img = img.compress(valid)
        # Precess the pixel coordinates to the (J2000) epoch of the map
        m = a.coord.convert_m('eq','eq', 
            iepoch=obs_date, oepoch=ephem.J2000)
        ex,ey,ez = n.dot(m, n.array([ex,ey,ez])) 
        # Put the data into the skymap
        skymap.add((ex,ey,ez), map_wgts, img)
        #save skymap
        skymap.to_fits(mapFn, clobber=True)

