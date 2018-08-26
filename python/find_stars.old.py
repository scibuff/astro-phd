#!/usr/bin/python3

# Find stars in a FITS image

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import skimage.feature as feature
from scipy import ndimage

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#if len(sys.argv) == 1:
#  sys.exit("Usage: find_stars.py infile.fits outfile.fits pixels.txt")
#  exit()
#elif len(sys.argv) == 4:
#  infits = sys.argv[1]
#  outfits = sys.argv[2]
#  pixfile = sys.argv[3]
#else:
#  sys.exit("Usage: find_stars.py infile.fits outfile.fits pixels.txt")
#  exit() 
#
 
#infits = '/mnt/hgfs/data/stars.fits'
#outfits = 'out.fits'
infits = '/mnt/hgfs/data/m1.fits.fz'
outfits = 'out.fits.fz'
outpng = 'out.png'
pixfile = 'pixels.txt'
header_index = 1
 
# Set true for displayed list
verboseflag = True  

# Set the image bias such that black is equal to bias
bias = 0.

# Set fractional limit (0 to 1) in list of extrema 
# Discard points below this limit
limcut = 0.5

# Open the fits file readonly by default and create an input hdulist
inlist = pyfits.open(infits) 

# Assign the input header 
inhdr = inlist[header_index].header

# Assign image data to numpy array 
img =  inlist[header_index].data.astype('float32')

# Filter the image so that blacker than black are not processed
img[img < bias] = bias

# Get the dimensions for later use
xsize, ysize = img.shape


# Filter the image to process only stars
lowpass = ndimage.gaussian_filter(img, 3)
highpass = img - lowpass

# Try an asinh transformation to handle a wide dynamic range
limg = np.arcsinh(highpass)
limg = limg / limg.max()
limg[limg < limcut] = bias

# Write this intermediate image so the user can see what is being done
outimage = limg
outlist = pyfits.PrimaryHDU(outimage.astype('float32'))
outlist.writeto(outfits, clobber = True)

# Find peaks in this filtered image
peaks = feature.peak_local_max(limg, min_distance=4, threshold_abs=limcut)

# Convert to FITS coordinates from numpy coordinates by adding 1
peaks = peaks + 1.

nx, ny = peaks.shape
nstars = nx

print("Found %d stars in %s" % (nstars, infits))

# Export FITS pixel coordinates to pixfile
pixfp = open(pixfile, 'a')

#Write the pixel table exchanging x and y 
#This provides the correct order for a FITS image
for i in range(nstars):
  xp, yp = peaks[i]
  pixline = "%5.3f  %5.3f \n" % (yp, xp)
  pixfp.write(pixline)

# Close the output file  
pixfp.close()

print('Mean:', np.mean(hdu.data))
print('Stdev:', np.std(hdu.data))
print('Min:', np.min(hdu.data))
print('Max:', np.max(hdu.data))


exit ()
