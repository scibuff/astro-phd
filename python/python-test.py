import platform
platform.python_version()

import astropy
import numpy as np
import photutils
import matplotlib

import os
result = os.system('clear')

print("astropy: " + astropy.__version__);
print("numpy: " + np.version.version);
print("photutils: " + photutils.version.version);
print("matplotlib: " + matplotlib.__version__);