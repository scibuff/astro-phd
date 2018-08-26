import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import aperture_photometry, CircularAperture, CircularAnnulus

def showStars( image ):
	
	stars = findStars(image)
	
	# take just the 5 brightest objects
	stars.sort('flux')
	stars.reverse()
	list = stars[:3]
	
	positions = zip(list['xcentroid'], list['ycentroid'])
	apertures = CircularAperture(positions, r=6) 
	annuli = CircularAnnulus(positions, r_in=18, r_out=24)
	phot_table = aperture_photometry(image, [apertures,annuli])
	
	plt.clf()
	
	plt.imshow( image, cmap='gray', norm=LogNorm() )
	annuli.plot(color='green', lw=1, alpha=0.75)
	apertures.plot(color='blue', lw=1, alpha=0.75)
	plt.savefig( '%s.png' % path, dpi=300 )
	
	print list
	print phot_table


##########################################################

import numpy as np
import astropy.io.fits as pyfits

from photutils import DAOStarFinder
from photutils import make_source_mask
from astropy.stats import median_absolute_deviation as mad
from astropy.stats import sigma_clipped_stats
from astropy.stats import mad_std
from math import sqrt
	
def findStars( image ):
	# In order to use the MAD as a consistent estimator for the estimation 
	# of the standard deviation σ, one takes σ = k * MAD where k is a constant 
	# scale factor, which depends on the distribution. For normally distributed
	# data k is taken to be k = 1.48[26]
	bkg_sigma = 1.48 * mad(image)
	t = 5*bkg_sigma
	daofind = DAOStarFinder(fwhm=3.0, threshold = t)
	stars = daofind(image)
	
	#stars['signal'] = stars['flux'] * t
	#
	data = image
	mask = make_source_mask(data, snr=2, npixels=5, dilate_size=11)
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
	madstd = mad_std(data);
	
	snrs = []
	for peak in stars['peak']:
		snrs.append( peak / madstd / 7.4 )
	stars['snr'] = snrs
	print((mean, median, std, bkg_sigma, t, madstd))
	#
	#print stars
	return stars

def processPath( path, index=0 ):
	
	image = getImage(path, index)
	showStars(image)

def getImage( path, index=0 ):
	
	fp = pyfits.open(path)
	hdu = fp[index]
	image = hdu.data.astype('float32')
	# H21 images are horizontally flipped
	image = np.fliplr(image) 
#	image -= np.median(image)
	return image

def getThreshold( path, index=0 ):
	fp = pyfits.open(path)
	hdu = fp[index]
	image = hdu.data.astype('float32')
	bkg_sigma = 1.48 * mad(image)
	return 5*bkg_sigma

images = [{'path': '/mnt/hgfs/data/k12h02h/00000338.2012FU62.REDUCED.FIT', 'dx': 0, 'dy': 0 },
	{'path': '/mnt/hgfs/data/k12h02h/00000339.2012FU62.REDUCED.FIT', 'dx': 0.03, 'dy': -0.67 },
	{'path': '/mnt/hgfs/data/k12h02h/00000340.2012FU62.REDUCED.FIT', 'dx': -0.05, 'dy': -0.59 },
	{'path': '/mnt/hgfs/data/k12h02h/00000341.CURSOR_POSITION.REDUCED.FIT', 'dx': -0.2, 'dy': -0.4 },
	{'path': '/mnt/hgfs/data/k12h02h/00000342.CURSOR_POSITION.REDUCED.FIT', 'dx': -0.07, 'dy': -0.55 },
	{'path': '/mnt/hgfs/data/k12h02h/00000343.CURSOR_POSITION.REDUCED.FIT', 'dx': -0.21, 'dy': -0.07 },
	{'path': '/mnt/hgfs/data/k12h02h/00000344.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.11, 'dy': -0.32 },
	{'path': '/mnt/hgfs/data/k12h02h/00000345.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.10, 'dy': -0.15 },
	{'path': '/mnt/hgfs/data/k12h02h/00000346.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.28, 'dy': -0.23 },
	{'path': '/mnt/hgfs/data/k12h02h/00000347.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.42, 'dy': -0.13 },
	#348 skipped
	{'path': '/mnt/hgfs/data/k12h02h/00000349.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.69, 'dy': 0.22 },
	{'path': '/mnt/hgfs/data/k12h02h/00000350.CURSOR_POSITION.REDUCED.FIT', 'dx': 0.85, 'dy': 0.55 },
	{'path': '/mnt/hgfs/data/k12h02h/00000351.CURSOR_POSITION.REDUCED.FIT', 'dx': 1.22, 'dy': 0.55 },
	{'path': '/mnt/hgfs/data/k12h02h/00000352.CURSOR_POSITION.REDUCED.FIT', 'dx': 1.32, 'dy': 1.00 },
	{'path': '/mnt/hgfs/data/k12h02h/00000353.CURSOR_POSITION.REDUCED.FIT', 'dx': 1.73, 'dy': 0.75 },
	{'path': '/mnt/hgfs/data/k12h02h/00000354.CURSOR_POSITION.REDUCED.FIT', 'dx': 1.81, 'dy': 0.96 }]

import scipy.ndimage as sciim

result = np.full((1024,1024),0,np.float32)

for row in images:
	
	path = row['path']
	index = 0
	
	fp = pyfits.open(path)
	hdu = fp[index]
	image = hdu.data.astype('float32')
	image = np.fliplr(image) 
	
	dx = row['dx']
	dy = row['dy']
	shifted = sciim.interpolation.shift( image, [dx, dy], order=1 )
	
	result += shifted

# flip the image back horizontally for H21
result = np.fliplr(result) 
result /= len(images)
	
outfile = '/mnt/hgfs/data/k12h02h/stack-%d.fits' % (len(images))
hdu = pyfits.PrimaryHDU( result )
hdu.writeto( outfile, overwrite=True )

##################
path = '/mnt/hgfs/data/k12h02h/00000339.2012FU62.REDUCED.FIT'
processPath(path)	
path = '/mnt/hgfs/data/k12h02h/stack-2.fits'
processPath(path)
path = '/mnt/hgfs/data/k12h02h/stack-4.fits'
processPath(path)
path = '/mnt/hgfs/data/k12h02h/stack-6.fits'
processPath(path)
path = '/mnt/hgfs/data/k12h02h/stack-9.fits'
processPath(path)
path = '/mnt/hgfs/data/k12h02h/stack-16.fits'
processPath(path)


##################
paths = ['/mnt/hgfs/data/k12h02h/00000347.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000349.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000350.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000351.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000352.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000353.CURSOR_POSITION.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/00000354.CURSOR_POSITION.REDUCED.FIT']
for path in paths:
	image = getImage(path)
	stars = findStars(image)
	stars.sort('flux')
	stars.reverse()
	list = stars[:3]
	print list
	print ''
	
##################
paths = ['/mnt/hgfs/data/k12h02h/00000339.2012FU62.REDUCED.FIT',
	'/mnt/hgfs/data/k12h02h/stack-2.fits',
	'/mnt/hgfs/data/k12h02h/stack-4.fits',
	'/mnt/hgfs/data/k12h02h/stack-6.fits',
	'/mnt/hgfs/data/k12h02h/stack-9.fits',
	'/mnt/hgfs/data/k12h02h/stack-16.fits',]
for path in paths:
	image = getImage(path)
	#print np.median(image), mad(image)
	bkg_sigma = 1.48 * mad(image)
	print mad(image), 5*bkg_sigma

##################
path = '/mnt/hgfs/data/k12h02h/00000339.2012FU62.REDUCED.FIT'
index = 0
image = getImage(path, index)
stars = findStars(image)
stars.sort('flux')
stars.reverse()
list = stars[:3]

positions = zip(list['xcentroid'], list['ycentroid'])
apertures = CircularAperture(positions, r=6) 
annuli = CircularAnnulus(positions, r_in=18, r_out=24)
phot_table = aperture_photometry(image, [apertures,annuli])