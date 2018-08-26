# https://reproject.readthedocs.io/en/stable/
# http://pythonhosted.org/pyplate/_modules/pyplate/solve.html
# http://123.physics.ucdavis.edu/observing_files/Python%20Image%20Reduction.pdf

import numpy as np
import scipy.ndimage as snd
import skimage.feature as feature
from astropy.io import fits

# Set the image bias such that black is equal to bias
bias = 0.0

#first = '/mnt/hgfs/data/stars.fits'

first = '/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0133-e91.fits.fz'
last = '/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0149-e91.fits.fz'

hdu = fits.open(first)
data = hdu[1].data

print('Mean:', np.mean(data))
print('Stdev:', np.std(data))

threshold = np.mean(data) + 5 * np.std(data)
print('Threshold:', threshold )

def find_max(data):
	xsize = len(data[0])
	ysize = len(data)
	maxx = 0
	maxy = 0
	max = data[maxx,maxy]
	for y in range(ysize):
		for x in range(xsize):
			if data[y,x] > max:
				maxx = x
				maxy = y
				max = data[y,x]
	return maxx, maxy

	
mask = data * 1
mask = np.multiply(mask,0.0)
y_size = len(data)
x_size = len(data[0])

for y in range(y_size):
	for x in range(x_size):
		if data[y,x] > threshold:
			mask[y,x] = 1

outfile = 'mask.fits'
hdu_out = fits.PrimaryHDU( mask )
hdu_out.writeto( outfile, clobber=True )

labels, num_features = snd.label( mask, np.ones((3,3)) )
centers = snd.center_of_mass( mask, labels, range( 1,  num_features + 1 ) )

x = np.array(centers)[:,0]
y = np.array(centers)[:,1]

coordinates = zip(y,x)
