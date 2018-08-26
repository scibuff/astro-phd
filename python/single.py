# http://docs.astropy.org/en/stable/wcs/

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle


#list = [ '/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0'+str(n)+'-e91.fits.fz' for n in range(133,149) ]
#concat = [ fits.getdata(image) for image in list ]
#final = np.sum(concat, axis=0)

fp = '/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0133-e91.fits.fz'
hdu = fits.open(fp)
hdu.info()
w = wcs.WCS(hdu[1].header)
print(w.wcs.name)
#w.wcs.print_contents()

pixcrd = np.array([[1024, 928]], np.float_)
world = w.wcs_pix2world(pixcrd, 1)

ra = Angle(world[0][0] * u.deg)
dec = Angle(world[0][1] * u.deg)
print( ra, dec )

# x and y are switched!c
#
# >>> f = fits.open( '/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0133-e91.fits.fz' )
# >>> d = f[1].data
# >>> print(d[1905,1731])   # get the pixel value at x=1732, y=1906
# 2100.68
# >>>

from urllib import urlencode
from urllib import urlretrieve


imsize = 20 #arcmin
ra_s = "%d:%d:%d" % (ra.hms.h, ra.hms.m, ra.hms.s )
dec_s = "%d:%d:%d" % (dec.dms.d, dec.dms.m, dec.dms.s )

# http://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=05:34:32&d=22:00:52&e=J2000&h=20&w=20&f=gif&c=none&fov=NONE
cutoutbaseurl = 'http://archive.stsci.edu/cgi-bin/dss_search'
query_string = urlencode(dict(f="gif", v="poss2ukstu_red", r=ra_s, d=dec_s, w=imsize, h=imsize, c="none", fov="NONE"))
url = cutoutbaseurl + '?' + query_string

# this downloads the image to your disk
urlretrieve(url, 'm1_cutout.jpg')