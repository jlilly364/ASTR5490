#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 18:31:53 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import time
from scipy.integrate import quad

# Class to model protoplanetary disk around Sun-like star
class Proplyd:
    
    def __init__(self,a=.1*10**-3,rho=2.0):
        # Inputs:
        #   a: mean grain radius (in m)
        #   rho: mean grain density (in g/cm^3 = kg/m^3)
        
        # Define global variables
        self.a = a*u.m
        self.rho = rho*u.kg/(u.m**3)
        
        # Calculate the dust sublimation radius using rad. equil. temp. eqn.
        self.r_sub = (const.R_sun/2*np.sqrt(1-0.3)*(5780/2000)**2).to(u.au)
        #print('r_sub = {0:.3f}'.format(r_sub.to(u.au)))
        
        # Calculate mean dust grain mass
        self.m_grain = 4/3*np.pi*self.a**3*self.rho
        #print(self.m_grain)
        print("{0} dust grains in this simulation".format(const.M_earth/self.m_grain))
        
Proplyd()
        