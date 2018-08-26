#!/usr/bin/python3

def addCircle( plt, ax, center ):
	y,x = center
	circle = plt.Circle( (x,y), 20, color='b', fill=False )
	ax.add_artist( circle )

def addStars( plt, ax, stars ):
	for star in stars:
		addCircle( plt, ax, star )

def findStars( hdu ):
	bias = 0.
	limcut = 0.8
	
	img = hdu.data.astype('float32')
	img[img < bias] = bias
	
	lowpass = ndimage.gaussian_filter(img, 3)
	highpass = img - lowpass
	
	limg = img;
	limg = np.arcsinh(highpass)
	limg = limg / limg.max()
	limg[limg < limcut] = bias
	
	peaks = feature.peak_local_max(limg, min_distance=4, threshold_abs=limcut)
	peaks = peaks + 1.
	return peaks

def processImage( path, index=0 ):
	fp = pyfits.open(path)
	hdu = fp[index]
	stars = findStars(hdu)
	
	plt.clf()
	fig, ax = plt.subplots()
	
	ax.imshow( hdu.data, cmap='gray', norm=LogNorm() )
	addStars( plt, ax, stars )
	
	plt.savefig( '%s.png' % path )

import numpy as np
import astropy.io.fits as pyfits
import skimage.feature as feature
from scipy import ndimage

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import (MinMaxInterval, SqrtStretch)
from astropy.visualization.mpl_normalize import ImageNormalize

#infits = '/mnt/hgfs/data/m1.fits.fz'
#infits = '/mnt/hgfs/data/stars.fits'
#outpng = 'stars.png'

processImage( '/mnt/hgfs/data/stars.fits', 0 )
processImage( '/mnt/hgfs/data/m1.fits.fz', 1 )
processImage( '/mnt/hgfs/data/ps1.fits', 0 )

exit ()
