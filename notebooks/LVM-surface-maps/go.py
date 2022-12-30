import numpy as np
import datetime
import glob

# Files created here need to be copied to LENS DESIGNS/From Others

def handle(filename):
    """ Read in a filename, write it out in a zemax sag format """
    dat = np.loadtxt(filename)
    x,y,z = dat.T

    xm = int(x.max()+1)
    ym = int(y.max()+1)

    dx = 280/538
    dy = 200/384
    z[~np.isfinite(z)] = 0
    z /= 1000 # micron -> mm

    #new = np.array(z.reshape(xm,ym)).T
    unitflag = 0 # Unit is mm (zemax manual)

    # From zemax manual, file format:
    # nx ny delx dely unitflag xdec ydec
    # Z dz/dx dz/dy d2z/dxdy nodata
    now = datetime.datetime.now()
    header = f"! Created on {now} using {__file__} \r\n{ym} {xm} {dx} {dy} {unitflag} 0 0"
    fmt = "%3.6f 0 0 0 0"
    print(header)
    np.savetxt(filename + ".dat", 
                z,
                fmt = fmt,
                header = header,
                comments = "",
                newline="\r\n")

def go():
    """ Find all the files and process them"""
    files = glob.glob("/Users/npk/Dropbox/REPOS/Konidaris-Research/notebooks/LVM-surface-maps/*.txt")
    print(files)

    for file in files:
        print(file)
        handle(file)


if __name__ == "__main__":
    go()