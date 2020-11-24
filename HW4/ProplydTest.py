#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 18:30:04 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u, constants as const
from MathTools import Planck
import time

wavelengths = np.linspace(10**-10,10**-4,10**4)*u.m
frequencies = const.c/wavelengths

def Loop():
    
    start_time = time.time()
    
    i = 0
    ydata = []
    for f in frequencies:
        i += 1
        helpme = Planck(f,5780*u.K)
        #print("Loop {0} of {1}: Planck at {2:.2e}={3:.2e}".format(i,len(frequencies),f,helpme))
        ydata.append(helpme.value)
    print("Loop took {0}".format(time.time()-start_time))
    """plt.plot(frequencies,ydata)
    plt.xscale('log')"""

def Call():
    start_time_new = time.time()
    
    test = Planck(frequencies,5780*u.K)
    print("Call took {0}".format(time.time()-start_time_new))
    #print(ydata-test)

Loop()
Call()