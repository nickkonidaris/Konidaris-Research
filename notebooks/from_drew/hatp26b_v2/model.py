from platon.transit_depth_calculator import TransitDepthCalculator
from platon.constants import M_jup, R_jup, R_sun, M_earth, R_earth

# All inputs and outputs for PLATON are in SI

# See https://arxiv.org/ftp/arxiv/papers/1705/1705.04354.pdf

Rs = 0.788 * R_sun
Mp = 18.6 * M_earth
Rp = 6.33 * R_earth

print((Rp / Rs)**2 / 1e-6)

# The initializer loads all data files.  Create a TransitDepthCalculator
# object and hold on to it
calculator = TransitDepthCalculator(method="xsec") #"ktables" for correlated k

# compute_depths is fast once data files are loaded
wavelengths, depths, info_dict = calculator.compute_depths(Rs, Mp, Rp, 903., logZ=1.96, CO_ratio=0.46, full_output=True, scattering_factor=10**2.37)
np.savetxt('model.dat', np.vstack([wavelengths, depths]).T)
