# Cold_sed_models

Following Carniani+19, the spectral energy distribution (SED) of the cold dust emission is modelled with a modified black-body (MBB) function given by:
 
$S_{\nu_{\rm obs}}^{\rm obs} = S_{\nu/(1+z)}^{\rm obs}$,  

and $S_{\nu/(1+z)}^{\rm obs} = \dfrac{\Omega}{(1+z)^3}(B_{\nu}(T_{\rm dust}(z))-B_{\nu}(T_{\rm CMB}(z)))(1-e^{-\tau_{\nu}})$, 

where $\Omega = (1+z)^4A_{\rm gal}D_{\rm L}^{-2}$ is the solid angle, $A_{\rm gal}$ and $D_{\rm L}$ are the surface area and luminosity distance of the galaxy, respectively. The dust optical depth is

$\tau_{\nu}=\dfrac{M_{\rm dust}}{A_{\rm galaxy}}k_0\biggl(\dfrac{\nu}{250\ \rm GHz}\biggr)^{\beta}$,

with $\beta$ the emissivity index, and $k_0 = 0.45\  \rm cm^{2}\ g^{-1}$ the mass absorption coefficient (Beelen+2006). The solid angle ($\Omega$) is estimated from the observed - or assumed - area of the galaxy. The effect of the CMB on the dust temperature is given by

$T_{\rm dust}(z)=((T_{\rm dust})^{4+\beta}+T_0^{4+\beta}[(1+z)^{4+\beta}-1])^{\frac{1}{4+\beta}}$,

with $T_0 = 2.73$ K.
We also considered the contribution of the CMB emission given by $B_{\nu}(T_{\rm CMB}(z)=T_0(1+z))$ (Da cunha+13).

# Packages needed

1. Astropy v 5.0.4
2. Numpy v 1.22.3
3. Scipy v 1.9.1
4. emcee v 3.1.1
5. matplotlib v 3.5.1
6. corner v 2.2.1
7. IPython v 7.22.0

# How to use the code

1. Download the repository in your working folder.
2. Open the Jupyter notebook named test.ipynb
3. In the notebook some examples to fit your data are reported
   
   3.1 There are three operational modes:
   - fit with 1 free parameter ($M_{\rm dust}$). Time for running the chain: ~ 1-2 min
   - fit with 2 free parameter ($M_{\rm dust}, \beta$). Time for running the chain: ~ 2-3 min
   - fit with 3 free parameter ($M_{\rm dust}, \beta, T_{\rm dust}$). Time for running the chain: ~ 3-4 min
     
   3.2 Along with some examples for the fitting procedure, one can find some examples for showing the best-fitting values (computed considering the 50th percentile of the distribution for each parameter) with their errors (i.e. 16th and 84th percentiles of the dirstribution for each parameter), and for plotting the chain, the corner plot, and the data with the best-fitting models.

   3.3 About the bayesian method to fit the data: we adoped a Gaussian likelihood with uniform priors. 4.0 < log(Mdust/Msun) < 9.0 and 0.5 < beta < 5.0 and 5. < Tdust/K < 300.

# Most important functions

- colddust_sed_models.sed_models.run_chain(start,npar,filename,nu_obs,flux,flux_err,z,Dl,A,fix,walker=0,iteration=3000,run=True,prog=True): run the chain
  
  - start = start position of the walkers in the chain (list of max 3 numbers)
    
  - npar = number of free parameter (it can be 1,2 or 3)
    
  - filename = filename used to save the chain (it has to be a string like 'path/name.h5'). Always remember to remove the file .h5 if you want to run twice the same chain.
    
  - nu_obs = observed frequencies of the data (in GHz)
    
  - flux = observed fluxes at nu_obs (in mJy)
    
  - flux_err = error associated to flux (in mJy)
    
  - z = redshift of the source
    
  - Dl = luminosity distance at z
    
  - A = area of the galaxy (in kpc^2)
    
  - fix = fixed parameters (It can be (beta, Tdust) or (Tdust) or none)
    
  - walker = number of walkers (Default: walker=0 that corresponds to walker = 10*npar). If one set walker=n=/0, then one has n walkers in the chain.
    
  - iteration = number of chains (Default: 3000. Usually enough for convergence).
    
  - run = condition to run the chain (Default: True). If set to False, the chain saved in filename will be imported.
    
  - prog = condition to show the progress of the MCMC (Default: True)
 
- colddust_sed_models.results_plot.sed_results(sampler, ndim, dis=250): display best-fitting values

   - sampler = chain made by run_chain() or saved in the .h5 file
     
   - ndim = number of free parameter (it can be 1,2 or 3); it is equal to npar in run_chain()
     
   - dis = number of chains to discard for burn-in (Default: 250)

- colddust_sed_models.results_plot.sed_res_plot(sampler, ndim, name, dis=250, path='', save=True, close=False): display and save the chain and the corner plot

  - sampler = as above
    
  - ndim = as above
    
  - name = string that will go into the filaname used to save the chain 'mcmc_chain_%s.png', and the corner plot 'corner_plot_%s.png'
    
  - dis = as before
    
  - path = path where to save the files (Default = '', i.e. working directory). If specified the path has to end with /, e.g. 'Users/Desktop/working/'
    
  - save = condition to save the plots in the working directory if path is not specified (Default = True).
    
  - close = condition to display the images in the notebook (Default = False, i.e. it displays the images)

 - colddust_sed_models.results_plot.plot_sed(freq_plot,sampler,nu_obs,flux,flux_err,z,Dl,A,c,fix,ndim,ylim,xlim,pos_text,name,dis=250,font='Times New Roman',path='',save=True,close=False): display and save the SED of the source with the best-fitting curve

   - freq_plot = x-axis array in observed frequency (GHz)
     
   - sampler = as above
     
   - nu_obs = as above
     
   - flux = as above
     
   - flux_err = as above
     
   - z, Dl, A = as above
     
   - c = velocity of light in m/s
     
   - fix = as above
     
   - ndim = as above
     
   - ylim = list of two elements (lower and upper limits for the y-axis)
     
   - xlim = list of two elements (lower and upper limits for the x-axis)
     
   - pos_text = list of two elements (x,y coordinates) for the position of the name of the source specified in name
     
   - name = string conteining the name of the source that will go in the pos_text and in the fileneme of savefig as 'sed-name-freq.png'
     
   - dis = as above
     
   - font = fontname for the axis and the text
     
   - path = as above
     
   - save = as above for saving the figure as 'sed-name-freq.png'
     
   - close = as above
    
