import os

def clear():
    n = os.system('clear')

clear()

import math
import numpy as np
import scipy.ndimage as sciim
import astropy.io.fits as pyfits
from astropy.table import Column
from astropy.stats import median_absolute_deviation as mad
from photutils import DAOStarFinder

def get_decimal_length():
    return 8

def get_rounded( f ):
    return round(f, get_decimal_length())

def get_rounded_format():
    return ".%df" % get_decimal_length()

# how many brightest stars to use from each image
def numberOfStars():
    return 15;

# returns a list of stars found by the DAOStarFinder
def getStars(image):
    # In order to use the MAD as a consistent estimator for the estimation
    # of the standard deviation σ, one takes σ = k * MAD where k is a constant
    # scale factor, which depends on the distribution. For normally distributed
    # data k is taken to be k = 1.48[26]
    bkg_sigma = 1.48 * mad(image)
    t = 5 * bkg_sigma
    daofind = DAOStarFinder(fwhm=4.0, threshold=t)
    stars = daofind(image)
    # print stars
    return stars


# loads a fits image given by the path
# index (optional) specified the hdu array index
def getImage(path, index=0):
    fp = pyfits.open(path)
    hdu = fp[index]
    image = hdu.data.astype('float32')
    # image -= np.median(image)
    return image

'''
# returns a list of stars for the purpose of being a list of reference stars, i.e. stars from the first image from
# a list of images for stacking
# the list is twice as long as a list of stars for regular images
def getRefStars(image):
    stars = getStars(image)
    list = stars[:numberOfStars() * 2]
    list.keep_columns(['id', 'xcentroid', 'ycentroid', 'flux'])
    return list


# returns the distance in pixels between star1 and star2
# the stars should be table rows returns from DAOStarFinder, i.e. from getStars or getRefStars routines
def getStarDistance(star1, star2):
    dx = star1['xcentroid'] - star2['xcentroid']
    dy = star1['ycentroid'] - star2['ycentroid']
    return math.sqrt(dx ** 2 + dy ** 2)


# takes a star, a list of reference stars and a limit radius (optional) and finds the
# reference star closest to the given star within that limit distance
# returns a tuple ( min distance, the closest reference star )
# the "closest reference star" will be None, if no ref. star is within the specified limit
def matchRefStar(star, ref, limit=10, offsetX=0, offsetY=0):
    dmin = 10 ** 6
    dminstar = None
    for refstar in ref:
        d = getStarDistance(star, refstar)
        if (d < dmin):
            dmin = d
            if (d < limit):
                dminstar = refstar
    return (dmin, dminstar)


# takes an astropy.Table from DAOStarFinder stars and matches them against a table of reference stars
# It adds columns called 'refid', 'd', 'dx', and 'dy', giving the reference star id, the total distances in px
# and the distance in x and y (of the centroids)
def matchRefStars(stars, ref, limit=10):
    refid = []
    dmins = []
    dx = []
    dy = []
    matched = stars.copy()
    for star in matched:
        # go thru the ref star list and find the closest match
        dmin, dminstar = matchRefStar(star, ref, limit)
        if (dminstar is not None):
            refid.append(dminstar['id'])
            dmins.append(dmin)
            dx.append(dminstar['xcentroid'] - star['xcentroid'])
            dy.append(dminstar['ycentroid'] - star['ycentroid'])
        else:
            refid.append(np.NaN)
            dmins.append(np.NaN)
            dx.append(np.NaN)
            dy.append(np.NaN)
    #
    matched['refid'] = Column(refid, description='Reference Star Id')
    matched['d'] = Column(dmins, description='Distance to the Reference Star')
    matched['dx'] = Column(dx, description='Distance to the Reference Star in x')
    matched['dy'] = Column(dy, description='Distance to the Reference Star in y')
    #
    return matched

hduIndex = 0
image = getImage(paths[0], hduIndex)
ref = getRefStars(image)
reg = []

result = np.full(image.shape, 0, np.float32)
avgx = 0
avgy = 0

for path in paths:
    image = getImage(path, hduIndex)
    stars = getStars(image)
    stars = stars[:numberOfStars()]
    stars.keep_columns(['id', 'xcentroid', 'ycentroid', 'flux'])
    # limit = 10 + math.sqrt( avgx**2 + avgy**2 )
    limit = 100
    matched = matchRefStars(stars, ref, limit)
    avgx = getDistanceAverage(matched['dx'])
    avgy = getDistanceAverage(matched['dy'])
    countNotNan = np.count_nonzero(~np.isnan(matched['dx']))
    reg.append([path, avgx, avgy, countNotNan])
    #
    # shifted = sciim.interpolation.shift(image, [-avgx, -avgy], order=1)
    # shifted = sciim.interpolation.shift(image, [0, 0], order=1)
    result += image

# flip the image back horizontally for H21
# result = np.fliplr(result)
result /= len(paths)

outfile = '/mnt/hgfs/data/k12h02h/stack-%d.fits' % (len(paths))
hdu = pyfits.PrimaryHDU(result)
hdu.writeto(outfile, overwrite=True)

print('path | avg(dx) | avg(dy) | count(non NaN)')
print('\n'.join('{}: {}'.format(*k) for k in enumerate(reg)))

# print(ref)
# print(stars)
# print(matched)
# print(reg)
######
'''

############################# angle functions ###########################
#   Returns the unit vector of the vector.
def unit_vector(vector):
    return vector / np.linalg.norm(vector)


#	Returns the angle in radians between vectors 'v1' and 'v2'::
#   >>> angle_between((1, 0, 0), (0, 1, 0))
#   1.5707963267948966
#   >>> angle_between((1, 0, 0), (1, 0, 0))
#   0.0
#   >>> angle_between((1, 0, 0), (-1, 0, 0))
#   3.141592653589793
def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


#	Returns a 2d vector from (x1,y1) to (x2,y2)
def make_vector(x1, y1, x2, y2):
    return [(x2 - x1), (y2 - y1)]


def make_star_vector(s1, s2):
    return make_vector(s1['xcentroid'], s1['ycentroid'], s2['xcentroid'], s2['ycentroid'])

########################## angle functions end ###########################

# takes a 'list' of values, converts it into an array via np.array() and returns the average
# of the values ignoring any np.NaN values
def getDistanceAverage(ds):
    a = np.array(ds)
    return np.average(np.ma.masked_array(a, np.isnan(a)))


# Takes a list of stars (from DAOStarFinder) and return a list of star vector angles. One star vector is the vector
# of the first two stars in the list (i.e. the two brightest stars) and the second star vector is the vector from the
# brightest star to the given star.
# Preconditions: the list of stars should be sorted based on the flux value in a descending order, i.e. the first
#                star in the list should be the brightest
def getStarAngles(stars):
    angles = [0, 0]
    star0 = stars[0]
    star1 = stars[1]
    v1 = make_star_vector(star0, star1)
    for i in range(2, len(stars), 1):
        star = stars[i]
        v2 = make_star_vector(star0, star)
        angle = angle_between(v1, v2)
        angles.append(angle)
    return angles

# Takes paths to images which will be stacked and add star angles column to the stars table.
#
def getStarsWithAngles( paths, hdu_index = 0 ):
    star_lists = []
    for path in paths:
        image = getImage(path, hdu_index)
        stars = getStars(image)
        columns = ['id', 'xcentroid', 'ycentroid', 'flux']
        for i in range( 1, len( columns ), 1 ):
            c = columns[i]
            stars[c].format = get_rounded_format()
        stars.keep_columns( columns )
        stars.sort('flux')
        stars.reverse()
        stars = stars[:numberOfStars()]
        angles = getStarAngles(stars)
        stars['angles'] = Column(angles, description='Angles')
        stars['angles'].format = get_rounded_format()
        #stars.sort('angles')
        star_lists.append(stars)
    return star_lists

# Takes a star and a list of reference stars and find the reference star with the closest reference angle
# return the angle diff, the reference star id, and the dx and dy distance between the star and the found reference
# star
def findStarWithClosestAngleMatch( star, refstars, limit = 0.0005 ):
    min = np.pi * 2 # max is 2 pi / or pi?
    dx, dy, min_star_id = ( np.NaN, np.NaN, np.NaN )
    for refstar in refstars:
        da = abs(refstar['angles'] - star['angles'])
        if ( da < limit ) and ( da < min ):
            min = da
            min_star_id = refstar['id']
            dx = refstar['xcentroid'] - star['xcentroid']
            dy = refstar['ycentroid'] - star['ycentroid']
    return min, min_star_id, dx, dy,

# Takes a star list and a reference star list and return 4 arrays containing:
# a) the angle difference between the star vectors from the star list and the reference star list
# b) the matched star id (from the reference list)
# c) pixel distance of the xcentroids (between a given star and the matched reference star)
# d) pixel distance of the ycentroids (between a given star and the matched reference star)
def getStarsWithReferenceMatch( stars, refstars ):
    das = [0, 0]
    matched_stars = [np.nan, np.nan]
    dxs = [0, 0]
    dys = [0, 0]
    for i in range(2, len(stars), 1):
        star = stars[i]
        da, matched_star_id, dx, dy = findStarWithClosestAngleMatch(star, refstars)
        das.append( get_rounded( da ) )
        matched_stars.append( matched_star_id )
        dxs.append( dx )
        dys.append( dy )
    return das, matched_stars, dxs, dys

# folder = '/mnt/hgfs/data/lcogtdata-20170112-196/'
# paths = [ folder + 'coj2m002-fs01-20161027-0133-e91.fits.fz',
#	folder + 'coj2m002-fs01-20161027-0134-e91.fits.fz',
#	folder + 'coj2m002-fs01-20161027-0135-e91.fits.fz']

folder = '/mnt/hgfs/data/k12h02h/'
paths = [folder + '00000338.2012FU62.REDUCED.FIT',
    folder + '00000364.CURSOR_POSITION.REDUCED.FIT',
    folder + '00000394.CURSOR_POSITION.REDUCED.FIT',
    folder + '00000420.CURSOR_POSITION.REDUCED.FIT']

hduIndex = 0 # H21 hdu index
#hduIndex = 1 # LGOCT hdu index
starLists = getStarsWithAngles( paths, hduIndex )
clear()

referenceStarList = starLists[0]

for i in range(1, len(starLists), 1):
    starList = starLists[i]
    das, matchedStars, dxs, dys = getStarsWithReferenceMatch( starList, referenceStarList )
    starList['da'] = Column(das, description='Angle Differences')
    starList['matches'] = Column(matchedStars, description='Matched Stars Id')
    starList['dx'] = Column(dxs, description='x distance to the matched star')
    starList['dy'] = Column(dys, description='y distance to the matched star')
    columns = ['da', 'dx', 'dy']
    for j in range(1, len(columns), 1):
        c = columns[j]
        starList[c].format = get_rounded_format()

print(starLists)

#### img reg ####

folder = '/mnt/hgfs/data/lcogtdata-20170112-196/'

hdu_index = 1
image0 = getImage(paths[0], hdu_index)
result = np.full(image.shape, 0, np.float32)

for path in paths:
    image = getImage(path, hdu_index)
    dy, dx = image_registration.cross_correlation_shifts(image0,image)
    shift = [ -get_rounded(dx), -get_rounded(dy) ]
    shifted = sciim.interpolation.shift(image, shift, order=1)
    result += shifted

# flip the image back horizontally for H21
# result = np.fliplr(result)
result /= len(paths)

outfile = folder + 'new-stack-%d.fits' % (len(paths))
hdu = pyfits.PrimaryHDU(result)
hdu.writeto(outfile, overwrite=True)

############################# debug functions #########################

from photutils import aperture_photometry, CircularAperture
import matplotlib.patches as patches
import matplotlib

matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm


def plotStars(path, image, stars):
    positions = zip(stars['xcentroid'], stars['ycentroid'])
    apertures = CircularAperture(positions, r=12)
    phot_table = aperture_photometry(image, apertures)
    plt.clf()
    plt.imshow(image, cmap='gray', norm=LogNorm())
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.savefig('%s.png' % path, dpi=300)
