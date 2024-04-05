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


def calib_err(nu_obs, flux, err):

    cal_err = np.array([])

    for freq in nu_obs:
        if freq>=80. and freq<=116.:
            cal_err = np.append(cal_err, 0.02)

        elif freq>=125. and freq<163.:
            cal_err = np.append(cal_err, 0.03)

        elif freq>=163. and freq<211.:
            cal_err = np.append(cal_err, 0.04)

        elif freq>=211. and freq<275.:
            cal_err = np.append(cal_err, 0.06)

        elif freq>=275. and freq<373.:
            cal_err = np.append(cal_err, 0.07)

        elif freq>=385. and freq<500.:
            cal_err = np.append(cal_err, 0.08)

        elif freq>=602. and freq<720.:
            cal_err = np.append(cal_err, 0.09)

        elif freq>=787. and freq<950.:
            cal_err = np.append(cal_err, 0.1)

        else:
            cal_err = np.append(cal_err, 0.0)
            
    flux_err = np.sqrt(err**2+(flux*cal_err)**2)
    
    return flux_err


def type_free_par(ndim, type_par):

    if ndim==1:
        if type_par==0:
            part = 'Md'
            
        elif type_par==1:
            part = 'beta'
            
        elif type_par==2:
            part = 'Td'

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')

    elif ndim==2:
        if type_par==0:
            part = np.array(['Md','beta'])
            
        elif type_par==1:
            part = np.array(['Md','Td'])
            
        elif type_par==2:
            part = np.array(['beta','Td'])

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')


    elif ndim==3:
        part = 0

    else:
        print('Error: ndim must be 1 or 2 or 3. Integer type.')


    return part


def fixing_order(ndim,type_par):

    if ndim==1:
        if type_par==0:
            print('Beta and Td must be fixed. Keep this order. E.g., fix_par=[1.6,50]')
            
        elif type_par==1:
            print('LogMd and Td must be fixed. Keep this order. E.g., fix_par=[8.,50]')
            
        elif type_par==2:
            print('LogMd and beta must be fixed. Keep this order. E.g., fix_par=[8.,1.6]')

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')

    elif ndim==2:
        if type_par==0:
            print('Td must be fixed. E.g., fix_par=[50]')
            
        elif type_par==1:
            print('Beta must be fixed. E.g., fix_par=[1.6]')
            
        elif type_par==2:
            print('LogMd must be fixed. E.g., fix_par=[8]')

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')


    elif ndim==3:
        print('No parameter to be fixed. Therefore, fix_par=0')

    else:
        print('Error: ndim must be 1 or 2 or 3. Integer type.')


    return

def st_order(ndim,type_par):

    if ndim==1:
        if type_par==0:
            print('LogMd free. Example for start position of the walkers, start_pos=[8.]')
            
        elif type_par==1:
            print('Beta free. Example for start position of the walkers, start_pos=[1.6]')
            
        elif type_par==2:
            print('Td free. Example for start position of the walkers, start_pos=[40.]')

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')

    elif ndim==2:
        if type_par==0:
            print('LogMd, beta free. Keep the order. \nExample for start position of the walkers, start_pos=[8.,1.6]')
            
        elif type_par==1:
            print('LogMd, Td free. Keep the order. \nExample for start position of the walkers, start_pos=[8.,50.]')
            
        elif type_par==2:
            print('Beta, Td free. Keep the order. \nExample for start position of the walkers, start_pos=[1.6, 50.]')

        else:
            print('Error: wrong parameter type. par must be 0, or 1 or 2. Integer type')


    elif ndim==3:
        print('LogMd, beta, Td free. Therefore, start_pos=[8.,1.6,50.]. Keep the order.')

    else:
        print('Error: ndim must be 1 or 2 or 3. Integer type.')


    return






