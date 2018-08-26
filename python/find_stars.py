#!/usr/bin/python3

def addCircle( plt, ax, center, image ):
	y,x = center
	print("%5.3f  %5.3f" % (x,y))
	xc, yc = centroid( image, x, y, 10 )
	print("%5.3f  %5.3f" % (xc,yc))
	circle = plt.Circle( (x,y), 20, color='b', fill=False )
	ax.add_artist( circle )

def addStars( plt, ax, stars, image ):
	print(stars)
	for star in stars:
		addCircle( plt, ax, star, image )

def findStars( data ):
	bias = 0.
	limcut = 0.6
	
	img = data.astype('float32')
	img[img < bias] = bias
	
	lowpass = ndimage.gaussian_filter(img, 3)
	highpass = img - lowpass
	
	limg = img;
	limg = np.arcsinh(highpass)
	limg = limg / limg.max()
	limg[limg < limcut] = 0.
	
	plt.clf()
	fig, ax = plt.subplots()
	ax.imshow( limg, cmap='gray')
	plt.savefig( 'findStars.png' )	
	
	peaks = feature.peak_local_max(limg, min_distance=4, threshold_abs=limcut)
	peaks = peaks + 1.
	return peaks

def processData( data, outpng ):	
	stars = findStars(data)
	plt.clf()
	fig, ax = plt.subplots()
	ax.imshow( data, cmap='gray', norm=LogNorm() )
	addStars( plt, ax, stars, data )
	plt.savefig( outpng )

def processImage( path, index=0 ):
	fp = pyfits.open(path)
	hdu = fp[index]
	processData(hdu.data,('%s.png' % path))


def cut(imdata, x, y, r):
	ysize, xsize = imdata.shape
	
	xmin = x-r
	xmax = x+r
	ymin = y-r
	ymax = y+r
	
	imin = int(np.floor(xmin - 0.5))
	imax = int(np.floor(xmax - 0.5))
	jmin = int(np.floor(ymin - 0.5))
	jmax = int(np.floor(ymax - 0.5))
	
	imin = max(0, imin)
	jmin = max(0, jmin)
	imax = min(xsize - 1, imax)
	jmax = min(ysize - 1, jmax)
	
	return imdata[jmin:jmax,imin:imax]
	
	

def centroid(imdata, x, y, r):
	if r < 0.5:
		return x, y
		
	ysize, xsize = imdata.shape
	
	xmin = x-r
	xmax = x+r
	ymin = y-r
	ymax = y+r
	
	imin = int(np.floor(xmin - 0.5))
	imax = int(np.floor(xmax - 0.5))
	jmin = int(np.floor(ymin - 0.5))
	jmax = int(np.floor(ymax - 0.5))
	
	imin = max(0, imin)
	jmin = max(0, jmin)
	imax = min(xsize - 1, imax)
	jmax = min(ysize - 1, jmax)
	
	tsig = 0.0
	txsig = 0.0
	tysig = 0.0
	r2 = r*r
	
	for i in range(imin, imax):
		for j in range(jmin, jmax):
			xp = float(i) + 1.0
			yp = float(j) + 1.0
			dx = xp - x
			dy = yp - y
			dx2 = dx * dx
			dy2 = dy * dy
			a2 = dx2 + dy2
			if a2 <= r2:
				txsig = dx * imdata[j, i] + txsig
				tysig = dy * imdata[j, i] + tysig
				tsig = imdata[j, i] + tsig
				
	tsig = max(tsig, 1.)
	xc = txsig/tsig + x
	yc = tysig/tsig + y
	
	return xc, yc

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
from mpl_toolkits.mplot3d import Axes3D

#infits = '/mnt/hgfs/data/m1.fits.fz'
#infits = '/mnt/hgfs/data/stars.fits'
#outpng = 'stars.png'

processImage( '/mnt/hgfs/data/stars.fits', 0 )
#processImage( '/mnt/hgfs/data/m1.fits.fz', 1 )
#processImage( '/mnt/hgfs/data/ps1.fits', 0 )

path = '/mnt/hgfs/data/stars.fits'
index = 0
fp = pyfits.open(path)
hdu = fp[index]
img = hdu.data
#181,314
#star = img[171:191,304:324]
star = img[151:211,284:344]
outpng = 'star.png'

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
x,y = np.meshgrid( np.arange(star.shape[1]), np.arange(star.shape[0] ) )
x = x.flatten()
y = y.flatten()
z = star.flatten()
ax.bar3d( x, y, np.zeros(len(z)), 1, 1, z )
proxy = plt.Rectangle((0,0),1,1,fc='b')
ax.legend([proxy],['peak [30,30]'])
plt.savefig(outpng)

exit ()
