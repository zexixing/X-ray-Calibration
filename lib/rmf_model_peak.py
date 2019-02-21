"""
define the model of instrumental response in astropy
"""
import os
import pylab as pl
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d
'''
def loss(x, l, f0_Ein):
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

def single_integ(m, n, A):
    '''
    integration of gaussian function
    '''
    new_A = A*np.sqrt(np.pi)/2.
    integ = erf(n) - erf(m)
    return new_A*integ

def single_integ_coef(E, Ein,
                      F, omiga, trap_n, read_n,
                      l, f0_Ein,
                      N_ph, N_flo, N_esc,
                      num_integ):
    '''
    integration coefficient x from 0 to l
    for flo(rescent), esc(ape peak) and lossed ph(oton)
    '''
    a = mass_attenuation(Ein)
    k = Ein*(1.-f0_Ein)/l
    b_ph = Ein*f0_Ein
    r_ph = resolution(Ein, F, omiga, trap_n, read_n)
    X_ph = np.square(a)*np.square(r_ph)/(2.*np.square(k))-a*(E-b_ph)/k
    A_ph = N_ph*(np.sqrt(2)*r_ph/k)*np.exp(X_ph)
    m_ph = a*r_ph/(np.sqrt(2)*k) - (E-b_ph)/(np.sqrt(2)*r_ph)
    n_ph = a*r_ph/(np.sqrt(2)*k) + (k*l-E+b_ph)/(np.sqrt(2)*r_ph)
    if X_ph >= 709.782 or erf(m_ph) - erf(n_ph) == 0.:
        num_integ = True
    if Ein > 1.74:
        b_esc = Ein*f0_Ein - 1.74
        r_esc = resolution(Ein-1.74, F, omiga, trap_n, read_n)
        X_esc = np.square(a)*np.square(r_esc)/(2.*np.square(k))-a*(E-b_esc)/k
        A_esc = N_esc*(np.sqrt(2)*r_esc/k)*np.exp(X_esc)
        m_esc = a*r_esc/(np.sqrt(2)*k) - (E-b_esc)/(np.sqrt(2)*r_esc)
        n_esc = a*r_esc/(np.sqrt(2)*k) + (k*l-E+b_esc)/(np.sqrt(2)*r_esc)
        if X_esc >= 709.782 or erf(m_esc) - erf(n_esc) == 0.:
            num_integ = True
        r_flo = resolution(1.74, F, omiga, trap_n, read_n)
        A_flo = N_flo\
                *np.exp(-0.5*np.square((E-1.74)/r_flo))\
                *(1.-np.exp(-a*l))/a
    else:
        A_esc, m_esc, n_esc, A_flo = 0., 0., 0., 0.
    integ_coef = {}
    integ_coef['A_ph'] = A_ph
    integ_coef['m_ph'] = m_ph
    integ_coef['n_ph'] = n_ph
    integ_coef['A_esc'] = A_esc
    integ_coef['m_esc'] = m_esc
    integ_coef['n_esc'] = n_esc
    integ_coef['A_flo'] = A_flo
    return integ_coef   

def integrate(f, a, b, N):
    '''
        self-defined numerical integration
    '''
    x = np.linspace(a+(b-a)/(2*N), b-(b-a)/(2*N), N)
    area = f(x[0])*(b-a)/N
    for i in x[1:]:
        fx = f(i)
        area = area + fx*(b-a)/N
    return area

def single_num_integ(E, Ein, N, k, b, l, r):
    def num_integ_func(x):
        a = mass_attenuation(Ein)
        num_integ = N*np.exp(-a*x-np.square((E-(k*x+b))/(np.sqrt(2)*r)))
        return num_integ
    return integrate(num_integ, 0., l, 2000.)

def gauss(x, center, sigma, A):
    '''
    gaussian function
    '''
    g = A*np.exp(-0.5*np.square((x-center)/sigma))
    return g

def response(E, Ein,
             F, omiga, trap_n, read_n,
             N_ph_Ein, N_flo_Ein, N_esc_Ein):
    '''
    response without charge-loss
    Ein should be in the unit of keV
    '''
    main = gauss(E, Ein, resolution(Ein, F, omiga, trap_n, read_n), N_ph_Ein)
    if Ein > 1.74:
        fluoresence = gauss(E, 1.74,
                            resolution(1.74 , F, omiga, trap_n, read_n),
                            N_flo_Ein)
        escape = gauss(E, Ein-1.74,
                       resolution(Ein-1.74, F, omiga, trap_n, read_n),
                       N_esc_Ein)
    else:
        fluoresence = gauss(E, 1.74,
                            resolution(1.74, F, omiga, trap_n, read_n),
                            0.)
        escape = gauss(E, Ein-1.74,
                       resolution(1.74, F, omiga, trap_n, read_n),
                       0.)
    res = main + fluoresence + escape
    return res

def response_loss(E, Ein,
                  F, omiga, trap_n, read_n,
                  l, f0_Ein,
                  N_ph_Ein, N_flo_Ein, N_esc_Ein,
                  t):
    '''
    response modified by the surface charge loss
    '''
    num_integ = False
    coef = single_integ_coef(E, Ein,
                             F, omiga, trap_n, read_n,
                             l, f0_Ein,
                             N_ph_Ein, N_flo_Ein, N_esc_Ein,
                             num_integ)
    m_ph, n_ph, A_ph = coef['m_ph'], coef['n_ph'], coef['A_ph']
    m_esc, n_esc, A_esc = coef['m_esc'], coef['n_esc'], coef['A_esc']
    A_flo = coef['A_flo']
    if num_integ == False:
        res_loss_1 = single_integ(m_ph, n_ph, A_ph) \
                     + single_integ(m_esc, n_esc, A_esc) \
                     + A_flo
    else:
        if Ein > 1.74:
            res_loss_1 = single_num_integ(E, Ein, N_ph_Ein,
                                          Ein*(1.-f0_Ein)/l, Ein*f0_Ein, l,
                                          resolution(Ein, F, omiga, trap_n, read_n))\
                         + single_num_integ(E, Ein, N_esc_Ein,
                                            Ein*(1.-f0_Ein)/l, Ein*f0_Ein - 1.74, l,
                                            resolution(Ein-1.74, F, omiga, trap_n, read_n))\
                         + A_flo
        else:
            res_loss_1 = single_num_integ(E, Ein, N_ph_Ein,
                                          Ein*(1.-f0_Ein)/l, Ein*f0_Ein, l,
                                          resolution(Ein, F, omiga, trap_n, read_n))                                                
    res_loss_2 = (1/mass_attenuation(Ein)) \
                 *(np.exp(-l*mass_attenuation(Ein))-np.exp(-t*mass_attenuation(Ein)))\
                 *response(E, Ein,
                           F, omiga, trap_n, read_n,
                           N_ph_Ein, N_flo_Ein, N_esc_Ein)
    res_loss = res_loss_1 + res_loss_2
    return res_loss

def shelf(E, Ein,
          F, omiga, trap_n, read_n, N_sh_Ein):
    '''
    shelf
    '''
    r = resolution(Ein, F, omiga, trap_n, read_n)
    sh = single_integ(-E/(np.sqrt(2)*r),
                     (Ein-E)/(np.sqrt(2)*r),
                     N_sh_Ein/(np.sqrt(2)*np.pi))
    return sh

def response_rmf(E, Ein,
                 read_n,
                 l, f0_Ein,
                 N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                 t):
    '''
    rmf: photopeak + fluoresence + escape + charge loss( surface + shelf)
    '''
    F = 0.15
    omiga = 3.65E-03
    trap_n = 0.
    res_rmf = response_loss(E, Ein,
                            F, omiga, trap_n, read_n,
                            l, f0_Ein,
                            N_ph_Ein, N_flo_Ein, N_esc_Ein, 
                            t) \
              + shelf(E, Ein,
                      F, omiga, trap_n, read_n, N_sh_Ein)
    return res_rmf

#------functions below are used to calculate-----
#------the integration over E of rmf response-----
def double_integ_coef_1(E1, E2, Ein,
                        F, omiga, trap_n, read_n,
                        l, f0_Ein,
                        N_ph, N_flo, N_esc, num_integ):
    '''
    integration coefficient: x from 0 to l, E from E1 to E2
    for flo(rescent), esc(ape peak) and lossed ph(oton)
    '''
    a = mass_attenuation(Ein)
    k = Ein*(1-f0_Ein)/l
    b_ph = Ein*f0_Ein
    r_ph = resolution(Ein, F, omiga, trap_n, read_n)
#    c_ph = -N_ph*r_ph*np.sqrt(np.pi/2.)/a
#    X_ph = a*b_ph/k+np.square(a*r_ph/(np.sqrt(2)*k))*erf(1./np.sqrt(2.))*(k*l/r_ph)
#    delta_c_ph= np.exp(np.float128(X_ph))*(np.exp(-a*E2/k)-np.exp(-a*E1/k))
#    C_ph = c_ph*delta_c_ph
#    C_ph = np.float64(C_ph)
    if np.square(a*r_ph/k)/2.+a*b_ph/k < 709.
    c_ph = -N_ph*(r_ph/a)*np.sqrt(np.pi/2.)*np.exp(np.square(a*r_ph/k)/2.+a*b_ph/k)
    c1_ph = np.exp(-a*E2/k)*erf((1./np.sqrt(2.))*(a*r_ph/k+(k*l+b_ph)/r_ph-E2/r_ph))\
            - np.exp(-a*E1/k)*erf((1./np.sqrt(2.))*(a*r_ph/k+(k*l+b_ph)/r_ph-E1/r_ph))
    c2_ph = np.exp(-a*E2/k)*erf((1./np.sqrt(2.))*(a*r_ph/k+b_ph/r_ph-E2/r_ph))\
            - np.exp(-a*E1/k)*erf((1./np.sqrt(2.))*(a*r_ph/k+b_ph/r_ph-E1/r_ph))
    C_ph = c_ph*(c1_ph-c2_ph)
    A1_ph = -N_ph*r_ph*np.exp(-a*l)*np.sqrt(2.)/a
    m1_ph = E1/(np.sqrt(2.)*r_ph)-(k*l+b_ph)/(np.sqrt(2.)*r_ph)
    n1_ph = E2/(np.sqrt(2.)*r_ph)-(k*l+b_ph)/(np.sqrt(2.)*r_ph)
    A2_ph = -N_ph*r_ph*np.sqrt(2.)/a
    m2_ph = E1/(np.sqrt(2.)*r_ph)-b_ph/(np.sqrt(2.)*r_ph)
    n2_ph = E2/(np.sqrt(2.)*r_ph)-b_ph/(np.sqrt(2.)*r_ph)
    if Ein > 1.74:
        b_esc = Ein*f0_Ein - 1.74
        r_esc = resolution(Ein-1.74, F, omiga, trap_n, read_n)
#        c_esc = -N_esc*r_esc*np.sqrt(np.pi/2.)/a
#        X_esc = a*b_esc/k+np.square(a*r_esc/(np.sqrt(2)*k))*erf(1./np.sqrt(2.))*(k*l/r_esc)
#        delta_c_esc= np.exp(np.float128(X_esc))*(np.exp(-a*E2/k)-np.exp(-a*E1/k))
#        C_esc = c_esc*delta_c_esc
#        C_esc = np.float64(C_esc)
        c_esc = -N_esc*(r_esc/a)*np.sqrt(np.pi/2.)*np.exp(np.float128(np.square(a*r_esc/k)/2.+a*b_esc/k))
        c1_esc = np.exp(-a*E2/k)*erf((1./np.sqrt(2.))*(a*r_esc/k+(k*l+b_esc)/r_esc-E2/r_esc))\
                - np.exp(-a*E1/k)*erf((1./np.sqrt(2.))*(a*r_esc/k+(k*l+b_esc)/r_esc-E1/r_esc))
        c2_esc = np.exp(-a*E2/k)*erf((1./np.sqrt(2.))*(a*r_esc/k+b_esc/r_esc-E2/r_esc))\
                - np.exp(-a*E1/k)*erf((1./np.sqrt(2.))*(a*r_esc/k+b_esc/r_esc-E1/r_esc))
        C_esc = c_esc*(c1_esc-c2_esc)
        A1_esc = -N_esc*r_esc*np.exp(-a*l)*np.sqrt(2.)/a
        m1_esc = E1/(np.sqrt(2.)*r_esc)-(k*l+b_esc)/(np.sqrt(2.)*r_esc)
        n1_esc = E2/(np.sqrt(2.)*r_esc)-(k*l+b_esc)/(np.sqrt(2.)*r_esc)
        A2_esc = -N_esc*r_esc*np.sqrt(2.)/a
        m2_esc = E1/(np.sqrt(2.)*r_esc)-b_esc/(np.sqrt(2.)*r_esc)
        n2_esc = E2/(np.sqrt(2.)*r_esc)-b_esc/(np.sqrt(2.)*r_esc)
        r_flo = resolution(1.74, F, omiga, trap_n, read_n)
        A_flo = N_flo*(1.-np.exp(-a*l))*np.sqrt(2)*r_flo/a
        m_flo = (E1-1.74)/(np.sqrt(2)*r_flo)
        n_flo = (E2-1.74)/(np.sqrt(2)*r_flo)
    else:
        A1_esc, m1_esc, n1_esc = 0., 0., 0.
        A2_esc, m2_esc, n2_esc = 0., 0., 0.
        C_esc = 0.
        A_flo, m_flo, n_flo = 0., 0., 0.
    integ_coef = {}
    integ_coef['A1_ph'] = A1_ph
    integ_coef['m1_ph'] = m1_ph
    integ_coef['n1_ph'] = n1_ph
    integ_coef['A2_ph'] = A2_ph
    integ_coef['m2_ph'] = m2_ph
    integ_coef['n2_ph'] = n2_ph
    integ_coef['C_ph'] = C_ph
    integ_coef['A1_esc'] = A1_esc
    integ_coef['m1_esc'] = m1_esc
    integ_coef['n1_esc'] = n1_esc
    integ_coef['A2_esc'] = A2_esc
    integ_coef['m2_esc'] = m2_esc
    integ_coef['n2_esc'] = n2_esc
    integ_coef['C_esc'] = C_esc
    integ_coef['A_flo'] = A_flo
    integ_coef['m_flo'] = m_flo
    integ_coef['n_flo'] = n_flo
    return integ_coef

def double_integ_coef_2(E1, E2, Ein,
                        F, omiga, trap_n, read_n,
                        l, f0_Ein,
                        N_ph, N_flo, N_esc,
                        t):
    '''
    integration coefficient: x from l to t, E from E1 to E2
    for flo(rescent), esc(ape peak) and lossed ph(oton)
    '''
    a = mass_attenuation(Ein)
    r_ph = resolution(Ein, F, omiga, trap_n, read_n)
    A_ph = N_ph*(np.exp(-a*l)-np.exp(-a*t))*np.sqrt(2)*r_ph/a
    t1_ph = (E1-Ein)/(np.sqrt(2)*r_ph)
    t2_ph = (E2-Ein)/(np.sqrt(2)*r_ph)
    if Ein > 1.74:
        r_flo = resolution(1.74, F, omiga, trap_n, read_n)
        A_flo = N_flo*(np.exp(-a*l)-np.exp(-a*t))*np.sqrt(2)*r_flo/a
        t1_flo = (E1-1.74)/(np.sqrt(2)*r_flo)
        t2_flo = (E2-1.74)/(np.sqrt(2)*r_flo)
        r_esc = resolution(Ein-1.74, F, omiga, trap_n, read_n)
        A_esc = N_esc*(np.exp(-a*l)-np.exp(-a*t))*np.sqrt(2)*r_esc/a
        t1_esc = (E1-(Ein-1.74))/(np.sqrt(2)*r_esc)
        t2_esc = (E2-(Ein-1.74))/(np.sqrt(2)*r_esc)
    else:
        A_flo, t1_flo, t2_flo = 0., 0., 0.
        A_esc, t1_esc, t2_esc = 0., 0., 0.
    integ_coef = {}
    integ_coef['A_ph'] = A_ph
    integ_coef['m_ph'] = t1_ph
    integ_coef['n_ph'] = t2_ph
    integ_coef['A_flo'] = A_flo
    integ_coef['m_flo'] = t1_flo
    integ_coef['n_flo'] = t2_flo
    integ_coef['A_esc'] = A_esc
    integ_coef['m_esc'] = t1_esc
    integ_coef['n_esc'] = t2_esc
    return integ_coef

def double_integ_coef_3(E1, E2, Ein,
                        F, omiga, trap_n, read_n,
                        l, f0_Ein,
                        N):
    '''
    integration coefficient: E from E1 to E2
    for shelf
    '''
    a = mass_attenuation(Ein)
    r = resolution(Ein, F, omiga, trap_n, read_n)
    A = N*r/(2.*np.sqrt(np.pi))
    m1 = -E2/(np.sqrt(2)*r)
    n1 = -E1/(np.sqrt(2)*r)
    m2 = (Ein-E2)/(np.sqrt(2)*r)
    n2 = (Ein-E1)/(np.sqrt(2)*r)
    integ_coef = {}
    integ_coef['A'] = A
    integ_coef['m1'] = m1
    integ_coef['n1'] = n1
    integ_coef['m2'] = m2
    integ_coef['n2'] = n2
    return integ_coef

def erf_integ(m, n, A):
    '''
    analytic expression of the integration of erf function
    '''
    def body(x):
        return x*erf(x) + np.exp(-np.square(x))/np.sqrt(np.pi)
    return A*(body(n) - body(m))

def double_integ(E1, E2, Ein,
                 read_n,
                 l, f0_Ein,
                 N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                 t):
    '''
    integration result for rmf response func over E
    using analytic expression
    '''
    F = 0.15
    omiga = 3.65E-03
    trap_n = 0
    coef1 = double_integ_coef_1(E1, E2, Ein,
                                F, omiga, trap_n, read_n,
                                l, f0_Ein,
                                N_ph_Ein, N_flo_Ein, N_esc_Ein)
    coef2 = double_integ_coef_2(E1, E2, Ein,
                                F, omiga, trap_n, read_n,
                                l, f0_Ein,
                                N_ph_Ein, N_flo_Ein, N_esc_Ein,
                                t)
    coef3 = double_integ_coef_3(E1, E2, Ein,
                                F, omiga, trap_n, read_n,
                                l, f0_Ein,
                                N_sh_Ein)
    part1 = (single_integ(coef1['m1_ph'], coef1['n1_ph'], coef1['A1_ph'])
             - single_integ(coef1['m2_ph'], coef1['n2_ph'], coef1['A2_ph'])
             + coef1['C_ph']) + \
            (single_integ(coef1['m1_esc'], coef1['n1_esc'], coef1['A1_esc'])
             - single_integ(coef1['m2_esc'], coef1['n2_esc'], coef1['A2_esc'])
             + coef1['C_esc']) + \
            (single_integ(coef1['m_flo'], coef1['n_flo'], coef1['A_flo']))
    part2 = single_integ(coef2['m_ph'], coef2['n_ph'], coef2['A_ph'])\
            + single_integ(coef2['m_flo'], coef2['n_flo'], coef2['A_flo'])\
            + single_integ(coef2['m_esc'], coef2['n_esc'], coef2['A_esc'])
    part3 = erf_integ(coef3['m2'], coef3['n2'], coef3['A']) \
            - erf_integ(coef3['m1'], coef3['n1'], coef3['A'])
    double_integration = part1 + part2 + part3
    return double_integration

'''
read_n = 26.0156542,
l = 0.00014
t = 0.028
f0_Ein = 0.5
N_ph_Ein = 587.490007
N_flo_Ein = 0.
N_esc_Ein = 100.
N_sh_Ein = 1.3073e-08
Ein = 2.
E = np.arange(0., 2.5, 0.01)

y2 = response_rmf(E, Ein,
                  read_n,
                  l, f0_Ein,
                  N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                  t)
pl.plot(E, y2, 'k')
pl.show()
'''
