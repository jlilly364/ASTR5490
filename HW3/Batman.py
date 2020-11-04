#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:34:25 2020

@author: jimmy
"""
import batman
import numpy as np
import matplotlib.pyplot as plt

def BatmanModel(t0=0.0,P=1.0,rad_pl=0.1,a=15.0,i=87.0,e=0.0,w=90.0,coeff=[0.5,0.1],plot=False):
    
    # Initialize the model
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = t0#76.6764                    #time of inferior conjunction
    params.per = P#13.17562                       #orbital period
    params.rp = rad_pl#0.01602                       #planet radius (in units of stellar radii)
    params.a = a#13.8                        #semi-major axis (in units of stellar radii)
    params.inc = i#88.19                      #orbital inclination (in degrees)
    params.ecc = e#0.0                       #eccentricity
    params.w = w#90.                        #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = coeff#[0.4899, 0.1809]     #limb darkening coefficients [u1, u2]
    
    t = np.linspace(-.25, .25, 1000)  #times at which to calculate light curve
    t_new = np.interp(t, (t.min(), t.max()), (0, 0.2))
    m = batman.TransitModel(params, t)    #initializes model
    
    flux = m.light_curve(params)                    #calculates light curve
    
    if plot==True:
        fig,ax = plt.subplots()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.plot(t_new,flux)

    return (t_new,flux)
#BatmanModel()
