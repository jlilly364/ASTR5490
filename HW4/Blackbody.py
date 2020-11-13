#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:15:10 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
from MathTools import TrapIntegrate
import time

# Class to generate and analze spectral energy distributions (SEDs)
class SED:
    
    def __init__(self,Teff,lambda_min,lambda_max,N):
        # Inputs:
        #   T_eff: effective temperature of blackbody (in K)
        #   lambda_min: minimum wavelength to calc. Planck function over (in m)
        #   lambda_max: minimum wavelength to calc. Planck function over (in m)
        #   Number of subintervals to integrate over    
    
        self.Teff = Teff*u.K
        self.lambda_min = lambda_min*u.m
        self.lambda_max = lambda_max*u.m
        self.N = N
    
    # Function to calculate main part Planck function at given wavelength
    def Planck(self,wavelength):
        # Inputs:
        #   wavelength: wavelength to calculate 
        # Returns:
        #   value of Planck function at that wavelength
        
        # Calculate 2hc^2 (prefactor in Planck's function)
        prefactor = 2*const.h*const.c**2
        
        # Calculate hc/kT (constant in exponential of Planck's function)
        exp_factor = const.h*const.c/(const.k_B*self.Teff)

        B = prefactor*wavelength**(-5)/(np.exp(exp_factor/wavelength)-1)

        return(B)
    
    # Function to plot spectral energy distribution of star
    def SEDStar(self,plot=False):
        # Inputs:
        #   plot: boolean to decide to make plot of SED
        # Returns:
        #   Plot of star's SED
        
        # Determine when program began running
        start_time = time.time()
        
        # Define list of wavelengths
        x = np.linspace(self.lambda_min,self.lambda_max,self.N)
        
        # Calculate Planck function at each wavelength
        y = self.Planck(x)
        
        # Calculate wavelength at which blackbody peaks
        lambda_peak = (2.89*(10**3)*u.micron*u.K)/self.Teff
        lambda_peak = lambda_peak.to(u.m).value
    
        # Find where blackbody peaks from my calculations
        peak_loc = np.argmax(y)
        lambda_max = x[peak_loc]
        
        # Decide whether to plot SED or not
        if plot == True:
            
            plt.figure(figsize=(8,4))
            
            # Plot data with x-axis on log scale
            plt.plot(x,y)
            plt.xscale('log')
            
            # Plot wavelength peak
            plt.vlines(lambda_peak,min(y).value,max(y).value,colors='red',\
                       linestyles='dashed',label=r"(Analytical) $\lambda_{peak}=%.2e m$"%lambda_peak)
            plt.vlines(lambda_max.value,min(y).value,max(y).value,colors='red',\
                       linestyles='dashed',label=r'(Numerical) $\lambda_{max}=%.2e m$'%lambda_max.value)
            plt.legend()
            
            # Axes labels and titles
            plt.xlabel('Wavelength (m)',fontsize=14)
            plt.ylabel(r'Spectral Radiance ($\frac{J}{m^3 sr}$)',fontsize=14)
            plt.title(r'Spectral Energy Distribution for $T_{eff}=$%dK Star'%self.Teff.value,fontsize=18)
            plt.tight_layout()
        
        # Tell user how long function took to run
        print('Program took %.2f sec to run'%(time.time()-start_time))

"""test = SED(5780,10**-10,10**-4,10**5)
#test.Planck(10**-8)
test.SEDStar(True)"""