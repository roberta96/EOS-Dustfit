import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.interpolate
from scipy.integrate import quad

import corner

from IPython.display import display, Math


from .sed_models import *

def labels(ndim,type_par):

    if ndim==1:
        if type_par==0:
            label = ["\log(M_{dust}/M_\odot)"] 
            
        elif type_par==1:
            label = ["beta"]
            
        elif type_par==2:
            label = ["T_{dust}"]

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')

    elif ndim==2:
        if type_par==0:
            label = ["\log(M_{dust}/M_\odot)","beta"] 
            
        elif type_par==1:
            label = ["\log(M_{dust}/M_\odot)","T_{dust}"] 
            
        elif type_par==2:
            label = ["beta","T_{dust}"] 

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')


    elif ndim==3:
        label = ["\log(M_{dust}/M_\odot)","beta","T_{dust}"]

    else:
        print('Error: ndim must be 1 or 2 or 3. Integer type.')


    return label


def sed_results(sampler, ndim, par, dis=250):

    flat_samples = sampler.get_chain(discard=dis, flat=True)

    q = np.array([])
    popt = np.zeros(ndim)
    j = 0
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.append(q, np.diff(mcmc))
        popt[i] = mcmc[1]
        j = j+2
    #print(q) 
    #mdust = 10**(popt[0])
    #mdusterr_m = np.log(10)*mdust*q[0]
    #mdusterr_p = np.log(10)*mdust*q[1]

    cond = (ndim==1 and par==0) or (ndim==2 and par==0) or (ndim==2 and par==1)
    if cond:
        mdust_chain = 10**flat_samples[:, 0]
        mdust = np.percentile(mdust_chain, [16, 50, 84])
        qdust = np.diff(mdust)

        mdust_val = np.array([mdust[1], qdust[0], qdust[1]])

    return (popt, q, mdust_val) if cond else (popt, q)


def sed_res_plot(sampler, ndim, name, par_type, par, dis=250, path='', save=True, close=False):

    fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    #print(len(samples))

    cond = (ndim==1 and par==0) or (ndim==2 and par==0) or (ndim==2 and par==1)
    if cond:
        popt, q, mdust = sed_results(sampler, ndim, par, dis=dis)

    else:
        popt, q = sed_results(sampler, ndim, par, dis=dis)
    
    if ndim==1:
        if par_type == 'Md':
            labels = [r"$\log(M_{\rm dust}/M_\odot)$"]

        elif par_type == 'beta':
            labels = [r"$\beta$"]
            
        elif par_type == 'Td':
            labels = [r"$T_{\rm dust}$"]
        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')
            
        truth = popt
        
    elif ndim==2:
        if par_type[0] == 'Md' and par_type[1] == 'beta':
            labels = [r"$\log(M_{\rm dust}/M_\odot)$", r"$\beta$"]
            
        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            labels = [r"$\log(M_{\rm dust}/M_\odot)$", r"$T_{\rm dust}$"]
            
        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            labels = [ r"$\beta$", r"$T_{\rm dust}$"]
            
        else:
            print('Error: wrong parameter types. Should be [Md, beta], or [Md, Td], or [beta, Td]. String format. Case sensitive.')
            
        truth = [popt[0],popt[1]]
        
    else:
        labels = [r"$\log(M_{\rm dust}/M_\odot)$", r"$\beta$", r"$T_{\rm dust}$"]
        truth = [popt[0],popt[1],popt[2]]
        
    for i in range(ndim):
        if ndim==1:
            ax = axes
        else:
            ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    if ndim==1:
        axes.set_xlabel("step number")
    else:
        axes[-1].set_xlabel("step number")

    if save==True:
        plt.savefig(path+'mcmc_chain_%s.png' %name, dpi=800, bbox_inches='tight', pad_inches=0.01)
        
    if close==True:
        plt.close('all')

    flat_samples = sampler.get_chain(discard=dis, flat=True)
    
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
        plt.savefig(path+'corner_plot_%s.png' %name, dpi=800, bbox_inches='tight', pad_inches=0.01)

    if close==True:
        plt.close('all')

    return 

def plot_sed(freq_plot,sampler,par_type,par,nu_obs,flux,flux_err,z,Dl,A,c,fix,ndim,ylim,xlim,pos_text,name,filename,radio=0,norm_rad=0,display_res=False,pos_res=[0.1,0.9],dis=250,font='Times New Roman',path='',save=True,close=False):

    wave_plot = c*1e6/(freq_plot*1e9) #micrometri

    yd, yu = ylim
    xd, xu = xlim
    
    fig, ax = plt.subplots(figsize=(11,9))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis = 'both', which = 'major', direction = 'in', length = 15, color = 'k', top = True, right = True, labelsize = 15, width = 1.2)
    ax.tick_params(axis = 'both', which = 'minor', direction = 'in', length = 5, top = True, right = True,width = 1.2)
    plt.ylim(yd, yu)
    ax.set_xlim(xu,xd)
    
    ax.errorbar(nu_obs,flux, yerr = flux_err, capsize = 3, lw = 1.4, marker = 'd', ms =6, ls='none',color = 'cyan', markeredgecolor='mediumblue',ecolor='mediumblue', zorder =3)

    cond = (ndim==1 and par==0) or (ndim==2 and par==0) or (ndim==2 and par==1)
    if cond:
        popt, q, mdust = sed_results(sampler, ndim, par, dis=dis)

    else:
        popt, q = sed_results(sampler, ndim, par, dis=dis)

    inst = fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if ndim ==1:

        if par_type == 'Md':
            inst.beta_fix = fix[0]
            inst.T = fix[1]
            ax.plot(freq_plot,inst.sobs_1par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$\log(M/M_{\odot})=%.2f$' % (popt[0], ),
                r'$\beta=%.2f$ (fixed)' % (fix[0], ),
                r'$T_{\rm dust}=%i $K (fixed)' % (fix[1], )))

        elif par_type == 'beta':
            inst.dmass = fix[0]
            inst.T = fix[1]
            ax.plot(freq_plot,inst.sobs_1par_beta(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$\beta=%.2f$' % (popt[0], ),
                r'$T_{\rm dust}=%i $K (fixed)' % (fix[0], ),
                r'$\log(M/M_{\odot})=%.2f$ (fixed)' % (fix[1], )))

        elif par_type == 'Td':
            inst.dmass = fix[0]
            inst.beta = fix[1]
            ax.plot(freq_plot,inst.sobs_1par_T(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$T_{\rm dust}=%i $K' % (popt[0], ),
                r'$\log(M/M_{\odot})=%.2f$ (fixed)' % (fix[0], ),
                r'$\beta=%.2f$ (fixed)' % (fix[1], )))
                

        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')
        
    elif ndim==2:
        if par_type[0] == 'Md' and par_type[1] == 'beta':
            inst.T = fix
            if radio == 0:
                ax.plot(freq_plot,inst.sobs_2par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)
            else:
                inst.norm = norm_rad
                ax.plot(freq_plot,inst.sobs_2par_radio(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$\log(M/M_{\odot})=%.2f$' % (popt[0], ),
                r'$\beta=%.2f $' % (popt[1], ),
                r'$T_{\rm dust}=%i $K (fixed)' % (fix[0], )))

        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            inst.beta_fix2 = fix
            ax.plot(freq_plot,inst.sobs_2par_Md_Td(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$\log(M/M_{\odot})=%.2f$' % (popt[0], ),
                r'$T_{\rm dust}=%i $K' % (popt[1], ),
                r'$\beta=%.2f$ (fixed)' % (fix[0], )))

        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            inst.dmass = fix
            ax.plot(freq_plot,inst.sobs_2par_beta_Td(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

            if display_res == True:
                textstr = '\n'.join((
                r'$\beta=%.2f$' % (popt[0], ),
                r'$T_{\rm dust}=%i $K' % (popt[1], ),
                r'$\log(M/M_{\odot})=%.2f$ (fixed)' % (fix[0], )))


    else:
        ax.plot(freq_plot,inst.sobs_3par(freq_plot, popt), label = 'MBB', c = 'mediumblue', lw = 1.5, zorder =1)

        if display_res == True:
                textstr = '\n'.join((
                r'$\log(M/M_{\odot})=%.2f$' % (popt[0], ),
                r'$\beta=%.2f $' % (popt[1], ),
                r'$T_{\rm dust}=%i $K' % (popt[2], )))

    props = dict(boxstyle='round', facecolor='white', edgecolor = 'k')

    px1, py1 = pos_res
    ax.text(px1, py1, textstr, transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)
    
    font_dic = {'family':font, 'size':'x-large'}    
    plt.xlabel('Obs Frequency (GHz)', fontsize = 20, fontname = font)
    plt.ylabel(r'$\rm F_{\nu}$ (mJy)', fontsize = 20, fontname = font)

    px, py = pos_text
    plt.text(px,py,name,fontdict=font_dic, bbox=dict(edgecolor='k', facecolor='none')) 
    fig.tight_layout()

    if save==True:
        plt.savefig(path+'sed-%s-freq.png' %filename, dpi=800,bbox_inches='tight', pad_inches=0.01)

    if close==True:
        plt.close('all')

    return

def sed_func_plot(freq_plot,nu_obs,par_type,popt,ndim,z,A,Dl,fix,radio=0,norm_radio=0,dis=250):

    inst = fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if ndim ==1:

        if par_type == 'Md':
            inst.beta_fix = fix[0]
            inst.T = fix[1]
            res = inst.sobs_1par(freq_plot, popt)

        elif par_type == 'beta':
            inst.dmass = fix[0]
            inst.T = fix[1]
            res = inst.sobs_1par_beta(freq_plot, popt)
            
        elif par_type == 'Td':
            inst.dmass = fix[0]
            inst.beta = fix[1]
            res = inst.sobs_1par_T(freq_plot, popt)

        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')

    elif ndim==2:
        if par_type[0] == 'Md' and par_type[1] == 'beta':
            inst.T = fix
            if radio ==0:
                res = inst.sobs_2par(freq_plot, popt)
            else:
                inst.norm = norm_radio
                res = inst.sobs_2par_radio(freq_plot, popt)
        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            inst.beta_fix2 = fix
            res = inst.sobs_2par_Md_Td(freq_plot, popt)

        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            inst.dmass = fix
            res = inst.sobs_2par_beta_Td(freq_plot, popt)

        else:
            print('Error: wrong parameter types. Should be [Md, beta], or [Md, Td], or [beta, Td]. String format. Case sensitive.')

    else:
        res = inst.sobs_3par(freq_plot, popt)
    

    return res

def l_sfr(flat_samples,ndim,par_type,z,c,Dl,A,fix,a1,b1):

    inst = fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if ndim==1:

        if par_type == 'Md':
            inst.beta_fix = fix[0]
            inst.T = fix[1]

            I3_chain = np.array([])
            for ind in range(len(flat_samples)):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_1par(xint,flat_samples[ind])
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)

        elif par_type == 'beta':
            inst.dmass = fix[0]
            inst.T = fix[1]

            I3_chain = np.array([])
            for ind in range(len(flat_samples)):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_1par_beta(xint,flat_samples[ind])
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)
            
        elif par_type == 'Td':
            inst.dmass = fix[0]
            inst.beta = fix[1]

            I3_chain = np.array([])
            for ind in range(len(flat_samples)):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_1par_T(xint,flat_samples[ind])
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)

        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')


    elif ndim==2:

        if par_type[0] == 'Md' and par_type[1] == 'beta':
            inst.T = fix

            I3_chain = np.array([])
            for ind in range(len(flat_samples[:,0])):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_2par(xint,(flat_samples[ind,0],flat_samples[ind,1]))
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)
            
        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            inst.beta_fix2 = fix

            I3_chain = np.array([])
            for ind in range(len(flat_samples[:,0])):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_2par_Md_Td(xint,(flat_samples[ind,0],flat_samples[ind,1]))
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)

        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            inst.dmass = fix

            I3_chain = np.array([])
            for ind in range(len(flat_samples[:,0])):
                xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
                func = inst.sobs_2par_beta_Td(xint,(flat_samples[ind,0],flat_samples[ind,1]))
                I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint)) #(mJy*GHz)

        else:
            print('Error: wrong parameter types. Should be [Md, beta], or [Md, Td], or [beta, Td]. String format. Case sensitive.')


    else:

        I3_chain = np.array([])
        for ind in range(len(flat_samples[:,0])):
            xint = np.logspace(np.log10(b1),np.log10(a1),num=100)
            func = inst.sobs_3par(xint,(flat_samples[ind,0],flat_samples[ind,1],flat_samples[ind,2]))
            I3_chain = np.append(I3_chain,scipy.integrate.simps(func,xint))

    
    dl = Dl * 3.08e21 #cm
    Lir3 = 4*np.pi*(dl**2)*I3_chain*1e-26*1e9 #(cm2)*(mJy*GHz)*1e-26*1e9 = erg/s
    L_sun = 3.826*1e33 #erg/s
    LIR_chain = Lir3/L_sun #Lsun
    SFR_chain = 1e-10*LIR_chain #M_sun/year

    LIR = np.percentile(LIR_chain, [16, 50, 84])

    SFR = np.percentile(SFR_chain, [16, 50, 84])

    return LIR, SFR

    

def lfir_sfr(sampler,ndim,par_type,z,c,Dl,A,fix,dis=250):

    flat_samples = sampler.get_chain(discard=dis, flat=True)
    
    a1 = (c*1e-9/40e-6)/(1+z) #GHz
    b1 = (c*1e-9/1000e-6)/(1+z)

    LFIR, SFR = l_sfr(flat_samples,ndim,par_type,z,c,Dl,A,fix,a1,b1)
    
    return LFIR, SFR


def ltir_sfr(sampler,ndim,par_type,z,c,Dl,A,fix,dis=250):

    flat_samples = sampler.get_chain(discard=dis, flat=True)
    
    a1 = (c*1e-9/8e-6)/(1+z) #GHz
    b1 = (c*1e-9/1000e-6)/(1+z)

    LTIR, SFR = l_sfr(flat_samples,ndim,par_type,z,c,Dl,A,fix,a1,b1)

    return LTIR, SFR


def l_SFR_bestfit(popt,ndim,par_type,z,c,Dl,A,fix,a1,b1):

    inst = fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if ndim==1:

        arg = np.array([popt])

        if par_type == 'Md':
            inst.beta_fix = fix[0]
            inst.T = fix[1]
            I3 = quad(inst.sobs_1par, b1,a1, args=arg) #(mJy*GHz)

        elif par_type == 'beta':
            inst.dmass = fix[0]
            inst.T = fix[1]
            I3 = quad(inst.sobs_1par_beta, b1,a1, args=arg) #(mJy*GHz)
            
        elif par_type == 'Td':
            inst.dmass = fix[0]
            inst.beta = fix[1]
            I3 = quad(inst.sobs_1par_T, b1,a1, args=arg) #(mJy*GHz)

        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')


    elif ndim==2:

        arg = np.array([popt[0],popt[1]])

        if par_type[0] == 'Md' and par_type[1] == 'beta':
            inst.T = fix
            I3 = quad(inst.sobs_2par, b1,a1, args=arg) #(mJy*GHz)
            
        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            inst.beta_fix2 = fix
            I3 = quad(inst.sobs_2par_Md_Td, b1,a1, args=arg) #(mJy*GHz)

        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            inst.dmass = fix
            I3 = quad(inst.sobs_2par_beta_Td, b1,a1, args=arg) #(mJy*GHz)

        else:
            print('Error: wrong parameter types. Should be [Md, beta], or [Md, Td], or [beta, Td]. String format. Case sensitive.')


    else:

        arg = np.array([popt[0],popt[1],popt[2]])
        I3 = quad(inst.sobs_3par, b1,a1, args=arg) #(mJy*GHz)

    
    dl = Dl * 3.08e21 #cm
    Lir3 = 4*np.pi*(dl**2)*I3[0]*1e-26*1e9 #(cm2)*(mJy*GHz)*1e-26*1e9 = erg/s
    L_sun = 3.826*1e33 #erg/s
    LIR = Lir3/L_sun #Lsun
    SFR = 1e-10*Lir3/L_sun #M_sun/year

    return LIR, SFR


def lfir_sfr_bestfit(popt,ndim,par_type,z,c,Dl,A,fix):
    a1 = (c*1e-9/40e-6)/(1+z) #GHz
    b1 = (c*1e-9/1000e-6)/(1+z)

    LTIR, SFR = l_SFR_bestfit(popt,ndim,par_type,z,c,Dl,A,fix,a1,b1)

    return LTIR, SFR

def ltir_sfr_bestfit(popt,ndim,par_type,z,c,Dl,A,fix):
    a1 = (c*1e-9/8e-6)/(1+z) #GHz
    b1 = (c*1e-9/1000e-6)/(1+z)

    LTIR, SFR = l_SFR_bestfit(popt,ndim,par_type,z,c,Dl,A,fix,a1,b1)

    return LTIR, SFR



