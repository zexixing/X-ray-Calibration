"""
define the model of instrumental response in astropy
"""
import os
import pylab as pl
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

def loss(x, f0_Ein, l, c, tau):
    """
    loss function which describs charge losses
    reference: O.Godet, A.P.Beardmore, AA, 2009
    """
    A = l*(1-f0_Ein)/(l+tau*c)
    B = tau*c*(1-f0_Ein)/(l+tau*c)
    if x < l:
        y = f0_Ein + A*np.power(x/l, c)
    else:
        y = 1 - B*np.exp((l-x)/tau)
    return y
'''
def loss(x, f0_Ein, l, c, tau):
    if x<l:
        y = x*(1-f0_Ein)/l + f0_Ein
    else:
        y = 1
    return y
'''
def resolution(Ein, F, omiga, trap_n, read_n):
    """
    trap_n is trapping noise
    read_n is reading noise
    reference: A.D.Short, NIMPR, 2002
    """
    sigma = omiga*np.sqrt(F*Ein/omiga
                          + np.square(trap_n)
                          + np.square(read_n)
                          )
    return sigma

def mass_attenuation(Ein):

    """
    calculate mu by interpolation from 0 to 15 keV
    data source: NIST
    """
    A = np.array([1.00000, 1.50000, 1.83890, 1.83891,
                  2.00000, 3.00000, 4.00000, 5.00000,
                  6.00000, 8.00000, 10.0000, 15.0000])
    B = np.array([1.570E+03, 5.355E+02, 3.092E+02, 3.192E+03,
                  2.777E+03, 9.784E+02, 4.529E+02, 2.450E+02,
                  1.470E+02, 6.468E+01, 3.389E+01, 1.034E+01])
    f = interp1d(A, B, fill_value='extrapolate')
    density = 2.330E+00
    return f(Ein)*density

def integrate(f, a, b, N):
    x = np.linspace(a+(b-a)/(2*N), b-(b-a)/(2*N), N)
    area = f(x[0])*(b-a)/N
    for i in x[1:]:
        fx = f(i)
        area = area + fx*(b-a)/N
    return area

def rmf_analytical(E, Ein, t,
                   f0_Ein, l, c, tau,
                   read_n, A,
                   ):
    """
    an analytical model of the spectral re-distribution
    E, Ein: keV, t: centimeter
    omiga: 3.65E-03keV, F: 0.15, readn: 5
    order of mag: t~200.0E-04, f0_Ein~0.6,
                  l~0.2E-04, c~2, tau~0.2E-04
    reference: A.D.Short, NIMPR, 2002
    """
    F = 0.15
    omiga = 3.65E-03
    trap_n = 0
    def pulse_height_position(x):
        y = (np.exp(-x*mass_attenuation(Ein))
             *np.exp(-np.square(E - Ein*loss(x, f0_Ein, l, c, tau))
                     /(2.0*np.square(resolution(Ein, F, omiga, trap_n, read_n)))
                     )
             )*A
        return y
    value = integrate(pulse_height_position, 0.0, t, 5000)
    return value
#    if isinstance(E, float) or isinstance(E, int):
#        integrated, _ = quad(pulse_height_position, a=0.0, b=t, args=(E))
#        return integrated
#    else:
#        def integrated(E):
#            integrated, _ = quad(pulse_height_position, a=0.0, b=t, args=(E))
#            return np.array(integrated)
#        def vectorizeInt():
#            integrateArray = []
#            for i in range(len(E)):
#                integrateArray.append(integrated(E[i]))
#            return np.array(integrateArray)
#        return vectorizeInt()

