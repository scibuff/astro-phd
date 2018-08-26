#!/usr/bin/python3

def findStars( image ):
	# In order to use the MAD as a consistent estimator for the estimation 
	# of the standard deviation σ, one takes σ = k * MAD where k is a constant 
	# scale factor, which depends on the distribution. For normally distributed
	# data k is taken to be k = 1.48[26]
	bkg_sigma = 1.48 * mad(image)
	stars = daofind(image, fwhm=3.0, threshold = 5*bkg_sigma)
	#print stars
	return stars

def processImage( path, index=0 ):
	fp = pyfits.open(path)
	hdu = fp[index]
	image = hdu.data.astype('float32')
	image -= np.median(image)
	stars = findStars(image)
#
	positions = zip(stars['xcentroid'], stars['ycentroid'])
	apertures = CircularAperture(positions, r=12) 
	phot_table = aperture_photometry(image, apertures)
	#print phot_table
#
	plt.clf()
	image = hdu.data.astype('float32')
	plt.imshow( image, cmap='gray', norm=LogNorm() )
	apertures.plot(color='blue', lw=1.5, alpha=0.5)	
	plt.savefig( '%s.png' % path, dpi=300 )

import numpy as np
import astropy.io.fits as pyfits
from photutils import datasets
from photutils import daofind
from astropy.stats import median_absolute_deviation as mad
from photutils import aperture_photometry, CircularAperture

import matplotlib.patches as patches
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

#infits = '/mnt/hgfs/data/m1.fits.fz'
#infits = '/mnt/hgfs/data/stars.fits'
#outpng = 'stars.png'

processImage( '/mnt/hgfs/data/stars.fits', 0 )
processImage( '/mnt/hgfs/data/m1.fits.fz', 1 )
processImage( '/mnt/hgfs/data/ps1.fits', 0 )


path = '/mnt/hgfs/data/m1.fits.fz'
index = 1
fp = pyfits.open(path)
hdu = fp[index]
image = hdu.data.astype('float32')
image -= np.median(image)

image = np.array([[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,10000,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0],
	[0,0,0,0,0,0,0,0,0,0]],np.int32)

import scipy.ndimage as sciim
shifted = sciim.interpolation.shift( image, [-0.01, -0.67], order=1 )
print shifted
	
exit ()
