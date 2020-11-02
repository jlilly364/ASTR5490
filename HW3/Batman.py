#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:34:25 2020

@author: jimmy
"""
from PeriodicityTools import Periodicity
import batman
import numpy as np
import matplotlib.pyplot as plt

def BatmanModel(t0,P,rad_pl,a,i,e,w,coeff):
    # Initialize the model
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = t0#0.0                    #time of inferior conjunction
    params.per = P#1.0                       #orbital period
    params.rp = rad_pl#0.1                       #planet radius (in units of stellar radii)
    params.a = a#13.8                        #semi-major axis (in units of stellar radii)
    params.inc = i#87.                      #orbital inclination (in degrees)
    params.ecc = e#0.0                       #eccentricity
    params.w = w#90.                        #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = coeff#[0.4899, 0.1809]     #limb darkening coefficients [u1, u2]
    
    
    t = np.linspace(-.025, .025, 1000)  #times at which to calculate light curve
    m = batman.TransitModel(params, t)    #initializes model
    
    flux = m.light_curve(params)                    #calculates light curve
    
    plt.plot(t,flux)