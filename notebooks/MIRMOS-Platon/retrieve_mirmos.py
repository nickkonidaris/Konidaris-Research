from __future__ import print_function

import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import corner

from platon.fit_info import FitInfo
from platon.retriever import CombinedRetriever
from platon.constants import R_sun, R_jup, M_jup

from platon.transit_depth_calculator import TransitDepthCalculator
from platon.constants import M_jup, R_sun, R_jup

# All quantities in SI
Rs = 1.16 * R_sun     #Radius of star
Mp = 1.00 * M_jup     #Mass of planet
Rp = 1.00 * R_jup      #Radius of planet
T = 1200.              #Temperature of isothermal part of the atmosphere


dat1 = np.loadtxt("atran.smo.12555.dat.txt")
dat2 = np.loadtxt("atran.smo.12560.dat.txt")
_, latm1, tatm1 = dat1.T
_, latm2, tatm2 = dat2.T

latm = np.concatenate([latm2, latm1])/1e6
tatm = np.concatenate([tatm2, tatm1])

# latm is the atmospheric wave at tatm (atmospehric throughput



#create a TransitDepthCalculator object and compute wavelength dependent transit depths
bins = []
l = 0.95
waves = [0.7]
bins = [[0.65, 0.75]]

while l < 2.45:
    dl = l/100
    bins.append([l, l+dl])
    waves.append(l)
    l += dl


BANDS = True

l = np.array(waves)
in_band =   ((l>0.5) & (l < .9)) | \
            ((l>0.972) & (l < 1.124)) | \
            ((l > 1.153) & (l < 1.352)) |\
            ((l>1.466) & (l<1.807)) |\
            ((l>1.921) & (l<2.404))
oob = l>0
bins = np.array(bins)/1e6
depth_calculator = TransitDepthCalculator()
if BANDS:
    roi = in_band
else:
    roi = np.ones_like(waves)==1

bins = bins[roi]

depth_calculator.change_wavelength_bins(bins)
wavelengths, transit_depths = depth_calculator.compute_depths(
    Rs, Mp, Rp, T, CO_ratio=0.2, cloudtop_pressure=1e6, add_scattering=True,
    scattering_slope=4, scattering_factor=10)

errors = np.zeros_like(transit_depths)
errors[:] = 30e-6

T_guess = 1200.


for ix,dl in enumerate(bins):
    roi = (latm>dl[0]) & (latm<=dl[1])
    trans = np.median(tatm[roi])
    print(trans)
    # if emis = 0 then ppm is = 20
    snr_max = 1/20e-6
    N_obs = snr_max**2 * trans

    errors[ix] = 1/np.sqrt(N_obs)

errors[0] = 20e-6

#create a Retriever object
retriever = CombinedRetriever()

#create a FitInfo object and set best guess parameters
fit_info = retriever.get_default_fit_info(
    Rs=Rs, Mp=Mp, Rp=Rp, T=T_guess,
    logZ=0, CO_ratio=0.53, log_cloudtop_P=4,
    log_scatt_factor=0, scatt_slope=4, error_multiple=1, T_star=6091.)

#Add fitting parameters - this specifies which parameters you want to fit
#e.g. since we have not included cloudtop_P, it will be fixed at the value specified in the constructor

fit_info.add_gaussian_fit_param('Rs', 0.02*R_sun)
fit_info.add_gaussian_fit_param('Mp', 0.04*M_jup)

# Here, emcee is initialized with walkers where R is between 0.9*R_guess and
# 1.1*R_guess.  However, the hard limit on R is from 0 to infinity.
fit_info.add_uniform_fit_param('Rp', 0.9*Rp, 1.1*Rp)

fit_info.add_uniform_fit_param('T', 300, 3000, 0.5*T_guess, 1.5*T_guess)
fit_info.add_uniform_fit_param("log_scatt_factor", 0, 5, 0, 1)
fit_info.add_uniform_fit_param("logZ", -1, 3)
fit_info.add_uniform_fit_param("log_cloudtop_P", -0.99, 5)
fit_info.add_uniform_fit_param("error_multiple", 0.5, 5)

#Use Nested Sampling to do the fitting
result = retriever.run_multinest(bins, transit_depths, errors, \
				None, None, None, \
                             fit_info, rad_method="xsec")
plt.savefig("best_fit.png")

np.save("chain.npy", result.chain)
np.save("logp.npy", result.lnprobability)

fig = corner.corner(result.flatchain,
                    range=[0.99] * result.flatchain.shape[1],
                    labels=fit_info.fit_param_names)
fig.savefig("emcee_corner.png")

