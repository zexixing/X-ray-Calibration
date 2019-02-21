import sys
sys.path.append("/anaconda2/lib/python2.7/site-packages/")
import _mypath
import os
import numpy as np
import pylab as pl
from lib.rmf_model_peak_app import double_integ, response_rmf
from lmfit import minimize, Parameters, report_fit
from scipy.interpolate import interp1d

thisdir = os.path.dirname(__file__)
N = {}
ratio = 6.082011581884892

def arf_S_E(arf_name, S_E_name):
    '''
    consider arf for S_E: ea*E_h
    save the multiplied results in a new file and return the file name
    '''
    S_E_to_thisdir_path = '../docs/casa'
    arf_to_thisdir_path = '../docs/crab' 
    S_E_file = S_E_name
    arf_file = arf_name
    S_E_path = os.path.join(thisdir, S_E_to_thisdir_path, S_E_file)
    arf_path = os.path.join(thisdir, arf_to_thisdir_path, arf_file)
    S_E = np.loadtxt(S_E_path)
    arf = np.loadtxt(arf_path)
    ener_S_E = S_E.T[0]
    ampl_S_E = S_E.T[1]
    flag_S_E = S_E.T[2]
    ener_arf = arf.T[0]/1000.
    area_arf = arf.T[1]
    def arf_interp(E):
        f = interp1d(ener_arf, area_arf, fill_value='extrapolate')
        return f(E)
    savename = 'arf_'+S_E_name
    save_path = os.path.join(thisdir, S_E_to_thisdir_path, savename)
    f = open(save_path, 'w')
    for i in range(0, len(ener_S_E)):
        arf_m_S_E = ampl_S_E[i]*arf_interp(ener_S_E[i])
        f.write(str(ener_S_E[i])+' '+str(arf_m_S_E)+' '+str(int(flag_S_E[i]))+'\n')
    f.close()
    return savename
    
def rmf_model(params, E, S_E):
    resolution = 1. #theoretical resolu when the spectra func is a Delata func
    d = int(len(np.array(S_E).flatten())/3.)
    read_n = params['read_n']
    l = params['l']
    t = params['t']
    if isinstance(E, float) or isinstance(E, int):
        pulse_height = np.zeros(1)
    else:
        pulse_height = np.zeros(len(E))
    if d == 1:
        delta_Ein = resolution
        flag = str(int(S_E[2]))
        f0_Ein = params['f0_Ein_'+flag]
        N_ph_Ein = params['N_ph_Ein_'+flag]
        N_flo_Ein = params['N_flo_Ein_'+flag]
        N_esc_Ein = params['N_esc_Ein_'+flag]
        N_sh_Ein = params['N_sh_Ein_'+flag]
        N_Ein = double_integ(0., 5., S_E[0],
                             read_n,
                             l, f0_Ein,
                             N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                             t)
        pulse_height = pulse_height + \
                       response_rmf(E, S_E[0],
                                    read_n,
                                    l, f0_Ein,
                                    N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein, 
                                    t)*S_E[1]*delta_Ein/N_Ein
#                                    t)*S_E[1]*delta_Ein*ratio/N_Ein
    else:
        delta_Ein = S_E[1][0] - S_E[0][0]
        for i in range(0, d):
            flag = str(int(S_E[i][2]))
            f0_Ein = params['f0_Ein_'+flag]
            N_ph_Ein = params['N_ph_Ein_'+flag]
            N_flo_Ein = params['N_flo_Ein_'+flag]
            N_esc_Ein = params['N_esc_Ein_'+flag]
            N_sh_Ein = params['N_sh_Ein_'+flag]
            N_Ein = double_integ(0., 5., S_E[i][0],
                                 read_n,
                                 l, f0_Ein,
                                 N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                                 t)
            N['N_Ein_'+flag+'_'+str(int(i))] = N_Ein
            pulse_height = pulse_height + \
                           response_rmf(E, S_E[i][0],
                                        read_n,
                                        l, f0_Ein,
                                        N_ph_Ein, N_flo_Ein, N_esc_Ein, N_sh_Ein,
                                        t)*S_E[i][1]*delta_Ein/N_Ein
#                                        t)*S_E[i][1]*delta_Ein*ratio/N_Ein
    return pulse_height

def rmf_error(params, E, E_h, S_E, eps_data):
    return (rmf_model(params, E, S_E) - E_h)/eps_data

def spec_fitting(E, E_h, S_E, eps_data):
    params = Parameters()
    params.add('read_n', value = 25., min = 0.0, max = 50.0) 
    params.add('l', value = 2.e-05, min = 0.0, max = 1.0E-03) 
    params.add('delta', value = 8.e-03, min = 0.0, max = 0.1) 
    params.add('t', expr = 'l + delta')
    d = int(len(np.array(S_E).flatten())/3.)
    if d == 1:
        flag = str(int(S_E[2]))
        params.add('f0_Ein_'+flag, value = 0.85, min = 0.0, max = 1.0)
        params.add('N_ph_Ein_'+flag, value = 500., min = 0.0, max = 5000000.)
        params.add('N_sh_Ein_'+flag, value = 5.E-03, min = 0.0, max = 10000.)
        if S_E[0]>1.839:
            params.add('N_esc_Ein_'+flag, value = 50., min = 0., max = 3000000.)
            params.add('N_flo_Ein_'+flag, value = 50., min = 0., max = 5000000.)
        else:
            params.add('N_esc_Ein_'+flag, value = 0., vary = False)
            params.add('N_flo_Ein_'+flag, value = 0., vary = False)
    else:
        for i in range(0, d):
            flag = str(int(S_E[i][2]))
            params.add('f0_Ein_'+flag, value = 0.85, min = 0.0, max = 1.0)
            params.add('N_ph_Ein_'+flag, value = 500., min = 0.0, max = 5000000.)
            params.add('N_sh_Ein_'+flag, value = 5.E-03, min = 0.0, max = 10000.)
            if S_E[i][0]>1.839:
                params.add('N_esc_Ein_'+flag, value = 50., min = 0., max = 3000000.)
                params.add('N_flo_Ein_'+flag, value = 50., min = 0., max = 5000000.)
            else:
                params.add('N_esc_Ein_'+flag, value = 0., vary = False)
                params.add('N_flo_Ein_'+flag, value = 0., vary = False)
    result = minimize(rmf_error, params, args=(E, E_h, S_E, eps_data),
                      method='ampgo')
    return result



data_to_thisdir_path = '../docs/casa/'
input_file = 'spec_input_all'
output_file = 'spec_output_all'
arf_name = 'ea.txt'

input_file = arf_S_E(arf_name, input_file)
input_path = os.path.join(thisdir, data_to_thisdir_path, input_file)
output_path = os.path.join(thisdir, data_to_thisdir_path, output_file)
S_E = np.loadtxt(input_path)
inp = np.loadtxt(output_path).T[0]
outp = np.loadtxt(output_path).T[1]
error = np.loadtxt(output_path).T[2]
result = spec_fitting(inp, outp, S_E, error)
params = dict(result.params.valuesdict())
print report_fit(result)
print N

x = np.arange(0., 5., 0.005)
#x = np.arange(min(inp), max(inp), 0.005)
delta = S_E[1][0] - S_E[0][0]
y = rmf_model(params, x, S_E)

fig = pl.figure()
ax = fig.add_subplot(1, 1, 1)
pl.plot(inp, outp, 'ko', MarkerSize=1)
ax.set_xscale('linear')
ax.set_yscale('linear')
pl.xlabel('keV')
pl.ylabel('counts s-1 kev-1')
#pl.xlim(min(inp), max(inp))
#pl.ylim(0.5*min(outp), 2.*max(outp))
pl.errorbar(inp, outp, yerr=[error,error],
            fmt='ko', MarkerSize=1, alpha = 0.3)
pl.plot(x, y, 'b-')
pl.title('fitting result of casA spectra (considering ARF)')
save_name = 'arf_specline_all.png'
save_path = os.path.join(thisdir, data_to_thisdir_path, save_name)
pl.savefig(save_path)
pl.show()
