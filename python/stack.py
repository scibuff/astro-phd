import math
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Column
from astropy.stats import median_absolute_deviation as mad
from photutils import DAOStarFinder

# how many brightest stars to use from each image
def numberOfStars(): 
	return 10;

# returns a list of stars found by the DAOStarFinder
def findStars( image ):
	# In order to use the MAD as a consistent estimator for the estimation 
	# of the standard deviation σ, one takes σ = k * MAD where k is a constant 
	# scale factor, which depends on the distribution. For normally distributed
	# data k is taken to be k = 1.48[26]
	bkg_sigma = 1.48 * mad(image)
	t = 5*bkg_sigma
	daofind = DAOStarFinder(fwhm=4.0, threshold = t)
	stars = daofind(image)
	stars.sort('flux')
	stars.reverse()
	#print stars
	return stars

# loads a fits image given by the path
# index (optional) specified the hdu array index
def getImage( path, index=0 ):
	fp = pyfits.open(path)
	hdu = fp[index]
	image = hdu.data.astype('float32')
	#image -= np.median(image)
	return image
	
def getRefStars( image ):
	stars = findStars( image )
	list = stars[:numberOfStars() * 2]
	list.keep_columns(['id','xcentroid','ycentroid','flux'])
	return list
	
#folder = '/mnt/hgfs/data/lcogtdata-20170112-196/'
#paths = [ folder + 'coj2m002-fs01-20161027-0133-e91.fits.fz',
#	folder + 'coj2m002-fs01-20161027-0134-e91.fits.fz',
#	folder + 'coj2m002-fs01-20161027-0135-e91.fits.fz']

folder = '/mnt/hgfs/data/k12h02h/'
paths = [ folder + '00000338.2012FU62.REDUCED.FIT',
	folder + '00000339.2012FU62.REDUCED.FIT',
	folder + '00000340.2012FU62.REDUCED.FIT']

index = 0
image = getImage( paths[0], index )
ref = getRefStars( image )
print(ref)

image = getImage( paths[1], index )
stars = findStars( image )
stars = stars[:numberOfStars()]
stars.keep_columns(['id','xcentroid','ycentroid','flux'])
print(stars)

def getStarDistance( star1, star2 ):
	dx = star1['xcentroid'] - star2['xcentroid']
	dy = star1['ycentroid'] - star2['ycentroid']
	return math.sqrt( dx**2 + dy**2 )

# takes a star, a list of reference stars and a limit radius (optional) and finds the
# reference star closes to the given star within that limit distance
# returns a tuple ( min distance, the closest reference star )
# the "closest reference star" will be None, if no ref. star is within the specified limit
def matchRefStar( star, ref, limit = 10 ):
	dmin = 10 ** 6
	dminstar = None
	for refstar in ref:
		d = getStarDistance(star, refstar)
		if (d < dmin):
			dmin = d
			if( d < limit ):
				dminstar = refstar
	return ( dmin, dminstar )

# takes an astropy.Table from DAOStarFinder stars and matches them against a table of reference stars
# It adds columns called 'refid', 'd', 'dx', and 'dy', giving the reference star id, the total distances in px
# and the distance in x and y (of the centroids)
def matchRefStars( stars, ref ):
	refid = []
	dmins = []
	dx = []
	dy = []
	matched = stars.copy()
	for star in matched:
		# go thru the ref star list and find the closest match
		dmin, dminstar = matchRefStar(star, ref)
		if ( dminstar is not None ):
			refid.append( dminstar['id'] )
			dmins.append( dmin )
			dx.append( dminstar['xcentroid'] - star['xcentroid'] )
			dy.append(dminstar['ycentroid'] - star['ycentroid'])
		else:
			refid.append(np.NaN)
			dmins.append(np.NaN)
			dx.append(np.NaN)
			dy.append(np.NaN)
	#
	matched['refid'] = Column( refid, description='Reference Star Id')
	matched['d'] = Column( dmins, description='Distance to the Reference Star')
	matched['dx'] = Column( dx, description='Distance to the Reference Star in x')
	matched['dy'] = Column( dy, description='Distance to the Reference Star in y')
	#
	return matched

# takes a 'list' of values, converts it into an array via np.array() and returns the average
# of the values ignoring any np.NaN values
def getDistanceAverage( ds ):
	a = np.array( ds )
	return np.average( np.ma.masked_array( a, np.isnan(a) ) )

matched = matchRefStars( stars, ref )
print(ref)
print(stars)
print(matched)

######

for path in paths:
	# lgoct has index = 1
	#index = 1
	# h21 has index = 0
	index = 0
	image = getImage(path,index)
	# H21 images are horizontally flipped
	image = np.fliplr(image)
	print( np.max(image) )
	stars = findStars(image)
	stars.sort('flux')
	stars.reverse()
	#list = stars[:numberOfStars()]
	list = stars[:7]
	print( list )
	plotStars( path, image, list )
	
	
############################# debug functions #########################

from photutils import aperture_photometry, CircularAperture
import matplotlib.patches as patches
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

def plotStars( path, image, stars ):
	positions = zip(stars['xcentroid'], stars['ycentroid'])
	apertures = CircularAperture(positions, r=12)
	phot_table = aperture_photometry(image, apertures)
	plt.clf()
	plt.imshow( image, cmap='gray', norm=LogNorm() )
	apertures.plot(color='blue', lw=1.5, alpha=0.5)
	plt.savefig( '%s.png' % path, dpi=300 )
