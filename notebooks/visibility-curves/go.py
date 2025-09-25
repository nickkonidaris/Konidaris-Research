import numpy as np
from pylab import *
import scipy as sp
from scipy.special import jv



def J1(x):
    return jv(1, x)

def V(rho, Theta_rad):
    """ Visibility of a star """
    return J1(rho * Theta_rad)/(rho * Theta_rad)


mas = 1/206265/1000.

def Vbin(rho, IA=1, IB=1, TA=1*mas, TB=1*mas, TD=1*mas):
    print("HERE")
    return IA * V(rho, TA) + IB * V(rho, TB) * exp(1j * rho * TD)/(IA+IB)

def Vbinsq(x, **rest):
    return np.abs(Vbin(x, **rest))**2

print("run")