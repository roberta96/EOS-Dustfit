import numpy as np
import matplotlib.pyplot as plt

from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy.visualization import quantity_support

from scipy.optimize import curve_fit

import scipy.interpolate
from scipy.integrate import quad

import emcee
import corner

#version 17/02/2026


class fitClass:
    def __init__(self):
        pass

    def sobs_1par(self,x,theta): #self.T, self.beta_fix, z,A
            p0 = theta #p0=esp

            M_dust = 1*10**p0 * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = self.T #K

            beta = self.beta_fix
            
            k_nu = 0.45*(x*(1+self.z)/250)**beta

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+beta)+T_0**(4+beta)*((1+self.z)**(4+beta)-1))**(1/(4+beta))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr


            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_1par_beta(self,x,theta): #self.T, self.dmass, z,A
            p0 = theta #p0=esp

            M_dust = 1*10**self.dmass * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = self.T #K

            beta = p0
            
            k_nu = 0.45*(x*(1+self.z)/250)**beta

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+beta)+T_0**(4+beta)*((1+self.z)**(4+beta)-1))**(1/(4+beta))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr


            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_1par_T(self,x,theta): #self.dmass, self.beta_fix, z,A
            p0 = theta #p0=esp

            M_dust = 1*10**self.dmass * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = p0 #K

            beta = self.beta_fix
            
            k_nu = 0.45*(x*(1+self.z)/250)**beta

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+beta)+T_0**(4+beta)*((1+self.z)**(4+beta)-1))**(1/(4+beta))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr


            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs
    
    def sobs_2par(self,x,theta): #self.T, z,A,D_l
            p0, p1 = theta #p0=esp, p1=beta

            M_dust = 1*10**p0 * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = self.T #K
            
            k_nu = 0.45*(x*(1+self.z)/250)**p1

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+p1)+T_0**(4+p1)*((1+self.z)**(4+p1)-1))**(1/(4+p1))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr
            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_3par(self,x,theta): #self.z,A
            p0, p1, p2 = theta #p0=esp, p1=beta, p2=T

            M_dust = 1*10**p0 * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = p2 #K
            
            k_nu = 0.45*(x*(1+self.z)/250)**p1

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+p1)+T_0**(4+p1)*((1+self.z)**(4+p1)-1))**(1/(4+p1))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr
            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_2par_Md_Td(self,x,theta): #self.beta, z,A,D_l
            p0, p1 = theta #p0=esp, p1=Td

            M_dust = 1*10**p0 * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)

            beta = self.beta_fix2
            
            k_nu = 0.45*(x*(1+self.z)/250)**beta

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (p1**(4+beta)+T_0**(4+beta)*((1+self.z)**(4+beta)-1))**(1/(4+beta))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr


            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_2par_beta_Td(self,x,theta): #self.dmass, z,A,D_l
            p0, p1 = theta #p0=esp, p1=Td

            M_dust = 1*10**self.dmass * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            k_nu = 0.45*(x*(1+self.z)/250)**p0

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (p1**(4+p0)+T_0**(4+p0)*((1+self.z)**(4+p0)-1))**(1/(4+p0))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr


            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy
            #print(sobs)  mJy
            
            return sobs

    def sobs_2par_radio(self,x,theta): #self.T, z,A,D_l
            p0, p1 = theta #p0=esp, p1=beta

            M_dust = 1*10**p0 * 1.99e30 *1e3
            sigma_dust = M_dust/(self.A*9e42)
            
            T_dust = self.T #K
            
            k_nu = 0.45*(x*(1+self.z)/250)**p1

            T_0 = 2.73
            T_cmb_z = T_0*(1+self.z)
            Bcmb = BlackBody(temperature=T_cmb_z*u.K)

            T_dust_z = (T_dust**(4+p1)+T_0**(4+p1)*((1+self.z)**(4+p1)-1))**(1/(4+p1))
            Bdust = BlackBody(temperature=T_dust_z*u.K)

            tau = sigma_dust*k_nu
            omega = (((1+self.z)**4)*self.A*self.D_l**(-2)) * u.sr
            s = omega*(Bdust(x*(1+self.z)*u.GHz)-Bcmb(x*(1+self.z)*u.GHz))*(1-np.exp(-tau))/(1+self.z)**3
            sobs = s.value * 1e26 #mJy

            f = 0.1
            a = -0.8
            x0 = 1.4/(1+self.z)
            radio = (1-f)*self.norm*(x/x0)**a + f*self.norm*(x/x0)**(-0.1)
            #print(sobs)  mJy
            
            return sobs+radio

def log_prior_1par(theta): #flat priors
    p0 = theta
    if 3.0 < p0 < 10.0:
        return 0.0
    return -np.inf

def log_prior_1par_beta(theta): #flat priors
    p0 = theta
    if 0.5 < p0 < 5.0:
        return 0.0
    return -np.inf

def log_prior_1par_T(theta): #flat priors
    p0 = theta
    if 5. < p0 < 300.0:
        return 0.0
    return -np.inf

def log_prior_2par(theta): #flat priors
    p0, p1 = theta
    if 3.0 < p0 < 10.0 and 0.5 < p1 < 5.0:
        return 0.0
    return -np.inf

def log_prior_3par(theta): #flat priors
    p0, p1, p2 = theta #p0=esp, p1=beta, p2=T
    if 3.0 < p0 < 10.0 and 0.5 < p1 < 5.0 and 5. < p2 < 300.:
        return 0.0
    return -np.inf

def log_prior_2par_Md_Td(theta): #flat priors
    p0, p1 = theta #p0=esp, p1=T
    if 3.0 < p0 < 10.0 and 5. < p1 < 300.:
        return 0.0
    return -np.inf

def log_prior_2par_beta_Td(theta): #flat priors
    p0, p1 = theta #p0=esp, p1=T
    if 0.5 < p0 < 5.0 and 5. < p1 < 300.:
        return 0.0
    return -np.inf

def log_likelihood_1par(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_1par(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_1par_beta(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_1par_beta(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_1par_T(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_1par_T(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_2par(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_2par(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_2par_radio(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_2par_radio(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_3par(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_3par(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_2par_Md_Td(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_2par_Md_Td(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_likelihood_2par_beta_Td(theta, nuobs, fl, fl_err,inst):
    err = fl_err

    fluxmod = inst.sobs_2par_beta_Td(nuobs, theta)
        
    logl = np.sum(np.log(1./np.sqrt(2*np.pi*err**2)) -(fl-fluxmod)**2/(2.*err**2))
    return logl

def log_posterior_1par(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_1par(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_1par(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_1par_beta(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_1par_beta(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_1par_beta(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_1par_T(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_1par_T(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_1par_T(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_2par(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_2par(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_2par(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_3par(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_3par(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_3par(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_2par_Md_Td(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_2par_Md_Td(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_2par_Md_Td(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_2par_beta_Td(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_2par_beta_Td(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_2par_beta_Td(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def log_posterior_2par_radio(theta, nuobs, fl, fl_err,inst):
    
    lp = log_prior_2par(theta)
        
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    
    logl = log_likelihood_2par_radio(theta, nuobs, fl, fl_err,inst)
    logpos = logl+lp
    
    if not np.isfinite(logpos):
        return -np.inf, -np.inf, -np.inf
    return logpos, lp, logl

def run_chain(start,npar,par_type,filename,nu_obs,flux,flux_err,z,Dl,A,fix,radio=0,norm_radio=0,walker=0,iteration=3000,run=True,prog=True):

    if walker == 0:
        walk = 10*npar
    else:
        walk = walker


    inst = fitClass()

    inst.z = z
    inst.D_l = Dl
    inst.A = A

    if npar==1:
        if par_type == 'Md':
            pos0 = start
            pos = [pos0] + 1e-2 * np.random.randn(walk, npar)
            log_posterior = log_posterior_1par
            inst.beta_fix = fix[0]
            inst.T = fix[1]

        elif par_type == 'beta':
            pos0 = start
            pos = [pos0] + 1e-2 * np.random.randn(walk, npar)
            log_posterior = log_posterior_1par_beta
            inst.dmass = fix[0]
            inst.T = fix[1]
            
        elif par_type == 'Td':
            pos0 = start
            pos = [pos0] + 1e-2 * np.random.randn(walk, npar)
            log_posterior = log_posterior_1par_T
            inst.dmass = fix[0]
            inst.beta_fix = fix[1]

        else:
            print('Error: wrong parameter type. Should be Md, or beta, or Td. String format. Case sensitive.')
            
    elif npar==2:
        if par_type[0] == 'Md' and par_type[1] == 'beta':
            pos0, pos1 = start
            pos = [pos0,pos1] + 1e-2 * np.random.randn(walk, npar)
            if radio==0:
                log_posterior = log_posterior_2par
            else:
                log_posterior = log_posterior_2par_radio
                inst.norm = norm_radio
            inst.T = fix[0]
        elif par_type[0] == 'Md' and par_type[1] == 'Td':
            pos0, pos1 = start
            pos = [pos0,pos1] + 1e-2 * np.random.randn(walk, npar)
            log_posterior = log_posterior_2par_Md_Td
            inst.beta_fix2 = fix[0]

        elif par_type[0] == 'beta' and par_type[1] == 'Td':
            pos0, pos1 = start
            pos = [pos0,pos1] + 1e-2 * np.random.randn(walk, npar)
            log_posterior = log_posterior_2par_beta_Td
            inst.dmass = fix[0]

        else:
            print('Error: wrong parameter types. Should be [Md, beta], or [Md, Td], or [beta, Td]. String format. Case sensitive.')

    else:
        pos0, pos1, pos2 = start
        pos = [pos0,pos1,pos2] + 1e-2 * np.random.randn(walk, npar)
        log_posterior = log_posterior_3par
        
    nwalkers, ndim = pos.shape

    # Set up the backend
    # Don't forget to clear it in case the file already exists

    run_the_chain = run

    if run_the_chain == True:
        backend = emcee.backends.HDFBackend(filename)
        backend.reset(nwalkers, ndim)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(nu_obs, flux, flux_err,inst), backend=backend)    
        sampler.run_mcmc(pos, iteration, progress=prog)
    else:
        sampler = emcee.backends.HDFBackend(filename)

    return sampler
    
