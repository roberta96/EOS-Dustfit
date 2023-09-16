import numpy as np
import matplotlib.pyplot as plt

import sed_models 
import results_plot 

z = 6.522
D_l = 64929500 #kpc
scale = 5.564 #kpc/"
c = 299792458 #m/s

A = np.array([0.350*0.327])*np.pi*scale**2/4 #kpc^2
A = np.mean(A)

nu_obs =  np.array([95.33111,245.6728]) #rest frame frequency, GHz
nu = nu_obs*(1+z)
wave_obs = c*1e6/(nu_obs*1e9) #micrometri

flux = np.array([0.094,2.034]) #mJy 
err1 = np.array([0.016,0.054]) #0.047 for the first data for gauss2, 0.016 for gauss
quad1 = np.array([0.02,0.07])
flux_err = np.sqrt(err1**2+(flux*quad1)**2)
Msun = 1.99e30 #kg

fix_par = 40 #K

start_pos = [8.,1.5]

ndim = 2 
sampler = sed_models.run_chain(start_pos, ndim,'j0224-chain-T40.h5',nu_obs,flux,flux_err,z,D_l,A,fix_par)
print('Chain done')

popt, q, mdust = results_plot.sed_results(sampler, ndim, dis=250)

label_new = ["\log(M_{dust}/M_\odot)", r"\beta"]
j=0
for i in range(ndim):
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{+{2:.3f}}}"
    txt2 = txt.format(popt[i], q[j], q[j+1], label_new[i])
    display(Math(txt2))
    j = j+2

print('M_dust/M_sun = %.2E - %.2E + %.2E' %(mdust[0], mdust[1], mdust[2]))
