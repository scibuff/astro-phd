import astropy
import numpy as np

from astropy.io import fits

filepath = "/mnt/hgfs/data/lcogtdata-20170112-196/coj2m002-fs01-20161027-0133-e91.fits.fz"
hdu_list = fits.open(filepath)
hdu_list.info()

image_data = hdu_list[1].data
print(type(image_data))
print(image_data.shape)

hdu_list.close()

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))
