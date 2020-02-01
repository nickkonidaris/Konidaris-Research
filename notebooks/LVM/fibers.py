import numpy as np


# Figure out how to allocate LVM fibers

def nfiber_in_hex(nring):
    " Number of fibers in a perfect hexagon with nring rings "
    return 1+3*nring*(nring+1)

nphot = 12 # number of standard star fibers
nspec =  3 # number of spectrographs


print("#rings #block #fib sky/block TOTAL #/spec")
for ngroups in np.arange(18,40):
    for fibers_per_group in np.arange(18,40):
        for nrings in np.arange(24., 27.):

            nscifiber = nfiber_in_hex(nrings) # Total number of science fibers in the IFU

            nfiber_per_spec = ngroups*fibers_per_group # Total number of fibers in the spectrogrpah
            nfibers = nfiber_per_spec * nspec # Total number of fibers in LVM
            nsky = nfibers - nphot - nscifiber # Total number of sky fibers
            nsky_per_group = nsky/nspec/ngroups # Number of sky fibers per group

            if nfibers > 2050: continue
            if nfibers < 1900: continue
            if nsky_per_group < 2: continue

            print("%6i %5i %4i       %4.2f %5i %6i " % (nrings, ngroups, fibers_per_group, nsky_per_group, nfibers, nfibers/3))



        
