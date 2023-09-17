# Cold_sed_models

Following Carniani+19, the spectral energy distribution (SED) of the cold dust emission is modelled with a modified black-body (MBB) function given by
 
$S_{\nu_{\rm obs}}^{\rm obs} = S_{\nu/(1+z)}^{\rm obs}$  and 

$S_{\nu/(1+z)}^{\rm obs} = \dfrac{\Omega}{(1+z)^3}(B_{\nu}(T_{\rm dust}(z))-B_{\nu}(T_{\rm CMB}(z)))(1-e^{-\tau_{\nu}})$, 

where $\Omega = (1+z)^4A_{\rm gal}D_{\rm L}^{-2}$ is the solid angle with $A_{\rm gal}$ , and $D_{\rm L}$ is the surface area and luminosity distance of the galaxy, respectively. The dust optical depth is

$\tau_{\nu}=\dfrac{M_{\rm dust}}{A_{\rm galaxy}}k_0\biggl(\dfrac{\nu}{250\ \rm GHz}\biggr)^{\beta}$,

with $\beta$ the emissivity index and $k_0 = 0.45\  \rm cm^{2}\ g^{-1}$ the mass absorption coefficient (Beelen+2006). The solid angle is estimated from the observed or assumed area of the galaxy. The effect of the CMB on the dust temperature is given by

$T_{\rm dust}(z)=((T_{\rm dust})^{4+\beta}+T_0^{4+\beta}[(1+z)^{4+\beta}-1])^{\frac{1}{4+\beta}}$,

with $T_0 = 2.73$ K.
We also considered the contribution of the CMB emission given by $B_{\nu}(T_{\rm CMB}(z)=T_0(1+z))$ (Da cunha+13).

# How to use the code

1. Download the repository in your working folder.
2. Open the Jupyter notebook named test.ipynb
3. In the notebook some examples to fit your data are reported
   
   3.1 There are three operational modes:
   - fit with 1 free parameter ($M_{\rm dust}$)
   - fit with 2 free parameter ($M_{\rm dust}, \beta$)
   - fit with 3 free parameter ($M_{\rm dust}, \beta, T_{\rm dust}$)
     
   3.2 Along with some examples for the fitting procedure, one can find some examples for showing the best-fitting values (computed considering the 50th percentile of the distribution for each parameter) with their errors (i.e. 16th and 84th percentiles of the dirstribution for each parameter), and for plotting the chain, the corner plot, and the data with the best-fitting models.
