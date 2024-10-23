# Most important functions

- cdsed_modelling.sed_models.run_chain(start,npar,filename,nu_obs,flux,flux_err,z,Dl,A,fix,walker=0,iteration=3000,run=True,prog=True): run the chain
  
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
 
- cdsed_modelling.results_plot.sed_results(sampler, ndim, dis=250): display best-fitting values

   - sampler = chain made by run_chain() or saved in the .h5 file
     
   - ndim = number of free parameter (it can be 1,2 or 3); it is equal to npar in run_chain()
     
   - dis = number of chains to discard for burn-in (Default: 250)

- cdsed_modelling.results_plot.sed_res_plot(sampler, ndim, name, dis=250, path='', save=True, close=False): display and save the chain and the corner plot

  - sampler = as above
    
  - ndim = as above
    
  - name = string that will go into the filaname used to save the chain 'mcmc_chain_%s.png', and the corner plot 'corner_plot_%s.png'
    
  - dis = as before
    
  - path = path where to save the files (Default = '', i.e. working directory). If specified the path has to end with /, e.g. 'Users/Desktop/working/'
    
  - save = condition to save the plots in the working directory if path is not specified (Default = True).
    
  - close = condition to display the images in the notebook (Default = False, i.e. it displays the images)

 - cdsed_modelling.results_plot.plot_sed(freq_plot,sampler,nu_obs,flux,flux_err,z,Dl,A,c,fix,ndim,ylim,xlim,pos_text,name,dis=250,font='Times New Roman',path='',save=True,close=False): display and save the SED of the source with the best-fitting curve

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
    
