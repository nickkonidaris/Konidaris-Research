
import numpy as np
from astropy.io import fits as FF
import os
import pylab as plt

files = os.listdir(".")


for fn in files:
    if not fn.endswith(".fits"): continue
    print(fn)

    dat = FF.open(fn)[0].data

    df = np.zeros_like(dat, dtype=np.float64)
    df[:] = np.nan

    ne0 = dat != 0
    df[ne0] = dat[ne0]

    spec = np.nanmedian(df, axis=1)

    np.savetxt(fn.rstrip(".fits") + ".txt", spec)


plt.plot(spec)

plt.show()