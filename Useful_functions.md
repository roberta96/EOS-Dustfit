# Most important functions

- cdsed_modelling.sed_models.run_chain(start,npar,par_type,filename,nu_obs,flux,flux_err,z,Dl,A,fix,radio=0,norm_radio=0,walker=0,iteration=3000,run=True,prog=True): run the chain
  
  - start = start position of the walkers in the chain (list of max 3 numbers)
    
  - npar = number of free parameter (it can be 1,2 or 3)
 
  - par_type = is a string or array of strings, depending on the value of npar. Case sensitive. If npar=1, then par_type=['Md'] or ['Td'] or ['beta']. If npar=2, then par_type=['Md', 'beta'] or ['Md', 'T'] or [beta, T]. If npar=3, you do not need this keyword, so set part_type=0.
    
  - filename = filename used to save the chain (it has to be a string like 'path/name.h5'). Always remember to remove the file .h5 if you want to run twice the same chain.
    
  - nu_obs = observed frequencies of the data (in GHz)
    
  - flux = observed fluxes at nu_obs (in mJy)
    
  - flux_err = error associated to flux (in mJy)
    
  - z = redshift of the source
    
  - Dl = luminosity distance at z
    
  - A = area of the galaxy (in kpc^2)
    
  - fix = fixed parameters (It can be (beta, Tdust) or (Tdust) or none)
 
  - radio = this is whether to include the radio part to the fitting, but it is still a working in progress update. So please do not change the default settings
 
  - norm_radio = same as above
    
  - walker = number of walkers (Default: walker=0 that corresponds to walker = 10*npar). If one set walker=n=/0, then one has n walkers in the chain.
    
  - iteration = number of chains (Default: 3000. Usually enough for convergence).
    
  - run = condition to run the chain (Default: True). If set to False, the chain saved in filename will be imported.
    
  - prog = condition to show the progress of the MCMC (Default: True)
 
- cdsed_modelling.results_plot.sed_results(sampler, ndim, par, dis=250): display best-fitting values

   - sampler = chain made by run_chain() or saved in the .h5 file
     
   - ndim = number of free parameter (it can be 1,2 or 3); it is equal to npar in run_chain()
 
   - par = selection of the fitted parameter (it can be 0,1,2). If ndim=1, sel_par=0 corresponds to logMdust; sel_par=1 to beta; sel_par=2 to Tdust. If ndim=2, par=0 corresponds to [logMdust, beta]; par=1 to [logMdust, T]; par=2 to [beta, T]. Same as sel_par in the test-eos.ipynb.
     
   - dis = number of chains to discard for burn-in (Default: 250)

- cdsed_modelling.results_plot.sed_res_plot(sampler, ndim, name, par_type, par, dis=250, path='', save=True, close=False): display and save the chain and the corner plot

  - sampler = as above
    
  - ndim = as above
    
  - name = string that will go into the filaname used to save the chain 'mcmc_chain_%s.png', and the corner plot 'corner_plot_%s.png'
 
  - par_type = corresponding free parameter/s (string format) to ndim and par
 
  - par = as above 
    
  - dis = as above
    
  - path = path where to save the files (Default = '', i.e. working directory). If specified the path has to end with /, e.g. 'Users/Desktop/working/'
    
  - save = condition to save the plots in the working directory if path is not specified (Default = True).
    
  - close = condition to display the images in the notebook (Default = False, i.e. it displays the images)

 - cdsed_modelling.results_plot.plot_sed(freq_plot,sampler,nu_obs,flux,flux_err,z,Dl,A,c,fix,ndim,ylim,xlim,pos_text,name,filename,dis=250,font='Times New Roman',path='',save=True,close=False): display and save the SED of the source with the best-fitting curve

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
     
   - name = string conteining the name of the source that will go in the pos_text
  
   - filename = string conteining the name used to save the image as 'sed-filename-freq.png'. You can put the path in this name or separated using the following keyword path
          
   - dis = as above
     
   - font = fontname for the axis and the text
     
   - path = as above
     
   - save = as above for saving the figure as 'sed-name-freq.png'
     
   - close = as above
    
- sed_func_plot(freq_plot,nu_obs,par_type,popt,ndim,z,A,Dl,fix,radio=0,norm_radio=0,dis=250): display a MBB function given some parameters
  
  - freq_plot, nu_obs, par_type = as above
 
  - popt = array with the best-fitting values for the fitted free paramaters. Or one can choose some arbitrary values. The physical quantities and the order in which one has to provide these values are set by ndim and par_type
 
  - z,A,Dl = as above
 
  - array with the value of the fixed parameter/s. The physical quantities and the order in which one has to provide these values are set by ndim and par_type
 
  - radio,norm_radio,dis = as above
