# <img src="/EOS-logo-final-circular.png" width="20%"> <br> EOS-Dustfit <br>

## Contents
1. [Description](#Description)
2. [Installation and setup](#Setup)
    - [Package requirements](#Package_need)
    - [How to install `EOS-Dustfit`](#Install)
    - [How to use the code](#Example)
4. [Modelling basics](#Basics)

## <a name="Description"></a>Description

`EOS-Dustfit` is a powerful tool for fitting the spectral energy distribution (SED) of the cold dust emission of galaxies in the millimeter regime, probed by ALMA and NOEMA telescope, and derive dust properties and star formation rate. Eos, as the goddess of dawn, symbolically link to the idea of cosmic or interstellar dust that interacts with light. Dust in space scatters and absorbs light, shaping how we observe distant stars and galaxies, much like how the atmosphere's dust scatters sunlight during dawn.

## <a name="Setup"></a>Installation and setup

### <a name="Package_need"></a>Package requirements

1. Astropy v 5.0.4
2. Numpy v 1.22.3
3. Scipy v 1.9.1
4. emcee v 3.1.1
5. matplotlib v 3.5.1
6. corner v 2.2.1
7. IPython v 7.22.0

### <a name="Install"></a>How to install `EOS-Dustfit`

1. Dowload the `colddust_sed_models` directory
2. Discover the path of the python package you are using (e.g. entering 'import sys; sys.executable' on the python terminal)
3. Open a terminal in the folder where colddust_sed_models has been dowloaded
4. Enter (remember to change `python-path`)
   ```
   python-path/python -m pip install -e colddust_sed_models
   ```
6. To verify that the package is correctly installed, download the test.ipynb notebook in a directory that is different from the one containing colddust_sed_models and then run the first cell of the test.ipynb notebook

### <a name="Example"></a>How to use the code

1. Install the package (see previous point)
2. Open the Jupyter notebook named test.ipynb
3. In the notebook some examples to fit your data are reported
   
   3.1 There are three operational modes:
   - fit with 1 free parameter. Time for running the chain: ~ 1-2 min
   - fit with 2 free parameter. Time for running the chain: ~ 2-3 min
   - fit with 3 free parameter. Time for running the chain: ~ 3-4 min
     
   3.2 Along with some examples for the fitting procedure, one can find some examples for showing the best-fitting values (computed considering the 50th percentile of the distribution for each parameter) with their errors (i.e. 16th and 84th percentiles of the dirstribution for each parameter), and for plotting the chain, the corner plot, and the data with the best-fitting models.

   3.3 About the bayesian method to fit the data: we adoped a Gaussian likelihood with uniform priors. 4.0 < log(Mdust/Msun) < 9.0 and 0.5 < beta < 5.0 and 5. < Tdust/K < 300.

## <a name="Basics"></a>Modelling of the dust emission

Following Carniani+19, the SED of the cold dust emission is modelled with a modified black-body (MBB) function given by:
 
$S_{\nu_{\rm obs}}^{\rm obs} = S_{\nu/(1+z)}^{\rm obs}$,  

and $S_{\nu/(1+z)}^{\rm obs} = \dfrac{\Omega}{(1+z)^3}(B_{\nu}(T_{\rm dust}(z))-B_{\nu}(T_{\rm CMB}(z)))(1-e^{-\tau_{\nu}})$, 

where $\Omega = (1+z)^4A_{\rm gal}D_{\rm L}^{-2}$ is the solid angle, $A_{\rm gal}$ and $D_{\rm L}$ are the surface area and luminosity distance of the galaxy, respectively. The dust optical depth is

$\tau_{\nu}=\dfrac{M_{\rm dust}}{A_{\rm galaxy}}k_0\biggl(\dfrac{\nu}{250\ \rm GHz}\biggr)^{\beta}$,

with $\beta$ the emissivity index, and $k_0 = 0.45\  \rm cm^{2}\ g^{-1}$ the mass absorption coefficient (Beelen+2006). The solid angle ($\Omega$) is estimated from the observed - or assumed - area of the galaxy. The effect of the CMB on the dust temperature is given by

$T_{\rm dust}(z)=((T_{\rm dust})^{4+\beta}+T_0^{4+\beta}[(1+z)^{4+\beta}-1])^{\frac{1}{4+\beta}}$,

with $T_0 = 2.73$ K.
We also considered the contribution of the CMB emission given by $B_{\nu}(T_{\rm CMB}(z)=T_0(1+z))$ (Da cunha+13).

