import numpy as np
import matplotlib.pyplot as plt

import scipy.interpolate
from scipy.integrate import quad

import corner

from IPython.display import display, Math


from .sed_models import *


def sed_results(sampler, ndim, dis=200):

    flat_samples = sampler.get_chain(discard=dis, flat=True)

    if ndim==1:
        label_new = ["\log(M_{dust}/M_\odot)"]
    elif ndim==2:
        label_new = ["\log(M_{dust}/M_\odot)", r"\beta"]
    else:
        label_new = ["\log(M_{dust}/M_\odot)", r"\beta", r"T_{dust}"]
    

    q = np.array([])
    popt = np.zeros(ndim)
    j = 0
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.append(q, np.diff(mcmc))
        popt[i] = mcmc[1]
        j = j+2
    #print(q) 
    mdust = 10**(popt[0])
    mdusterr_m = np.log(10)*mdust*q[0]
    mdusterr_p = np.log(10)*mdust*q[1]

    mdust_val = np.array([mdust, mdusterr_m, mdusterr_p])

    return popt, q, mdust_val 


def sed_res_plot(sampler, ndim, name, dis=200, save=True):

    fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    #print(len(samples))

    popt, mdust = sed_results(sampler, ndim, dis)
    
    if ndim==1:
        labels = [r"$\log(M_{\rm dust}/M_\odot)$"]
        truth = [popt]
    elif ndim==2:
        labels = [r"$\log(M_{\rm dust}/M_\odot)$", r"$\beta$"]
        truth = [popt[0],popt[1]]
    else:
        labels = [r"$\log(M_{\rm dust}/M_\odot)$", r"$\beta$", r"$T_{\rm dust}$"]
        truth = [popt[0],popt[1],popt[2]]
        
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number");

    if save==True:
        plt.savefig('mcmc_chain_%s.png' %name, dpi=800)

    fig.close()
    
    label_dic = {'size':15}

    fig = corner.corner(flat_samples, labels=labels, truths=truth, 
                        color = 'darkblue', truth_color = 'deepskyblue', label_kwargs = label_dic)

    axes = np.array(fig.axes).reshape((ndim, ndim))

    j = 0
    for i in range(ndim):
        #print(i)
        #print(j)
        ax = axes[i, i]
        ax.axvline(popt[i]-q[j], color="deepskyblue", ls = '--')
        ax.axvline(popt[i]+q[j+1], color="deepskyblue", ls = '--')
        j = j+2
        ax.minorticks_on()
        ax.tick_params(axis = 'both', which = 'major', direction = 'in', length = 6, 
                   color = 'k', labelsize = 12, width = 1.2)
        ax.tick_params(axis = 'both', which = 'minor', direction = 'in', length = 3,width = 1.2)

    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.minorticks_on()
            ax.tick_params(axis = 'both', which = 'major', direction = 'in', length = 6, 
                       color = 'k', top = True, right = True, labelsize = 12, width = 1.2)
            ax.tick_params(axis = 'both', which = 'minor', direction = 'in', length = 3,width = 1.2,
                          top = True, right = True)

    if save==True:
        plt.savefig('corner_plot_%s.png' %name, dpi=800)

    fig.close()

def plot_sed(freq_plot,sampler,nu_obs,flux,flux_err,z,Dl,A,fix,ndim,ylim,pos_text,name,dis=200,save=True):

    wave_plot = c*1e6/(freq_plot*1e9) #micrometri

    yd, yu = ylim
    
    fig, ax = plt.subplots(figsize=(7,7))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis = 'both', which = 'major', direction = 'in', length = 15, color = 'k', top = True, right = True, labelsize = 15, width = 1.2)
    ax.tick_params(axis = 'both', which = 'minor', direction = 'in', length = 5, top = True, right = True,width = 1.2)
    plt.ylim(yd, yi)
    
    plt.errorbar(nu_obs,flux, yerr = flux_err, capsize = 3, lw = 1.4, marker = 'd', ms =6, ls='none',color = 'cyan', markeredgecolor='mediumblue',ecolor='mediumblue', zorder =3)

    popt = sed_results(sampler, ndim, dis)

    inst = sed_models.fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if ndim ==1:
        inst.beta_fix = fix[0]
        inst.T = fix[1]
        plt.plot(freq_plot,inst.sobs_1par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)
    elif ndim==2:
        inst.T = fix
        plt.plot(freq_plot,inst.sobs_2par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)
    else:
        plt.plot(freq_plot,inst.sobs_3par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

    font_dic = {'family':font, 'size':'x-large'}    
    plt.xlabel('Obs Frequency (GHz)', fontsize = 20, fontname = font)
    plt.ylabel(r'$\rm F_{\nu}$ (mJy)', fontsize = 20, fontname = font)

    px, py = pos_text
    plt.text(px,py,name,fontdict=font_dic, bbox=dict(edgecolor='k', facecolor='none')) 
    fig.tight_layout()

    if save==True:
        plt.savefig('sed-%s-freq.png' %name, dpi=800)
