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
from scipy.integrate import quad
import time

# Class to generate and analze spectral energy distributions (SEDs)
class SED:
    
    def __init__(self,Teff=5780,lambda_min=10**-10,lambda_max=10**-4,N=10**5):
        # Inputs:
        #   T_eff: effective temperature of blackbody (in K)
        #   lambda_min: minimum wavelength to calc. Planck function over (in m)
        #   lambda_max: minimum wavelength to calc. Planck function over (in m)
        #   Number of subintervals to integrate over    
    
        self.Teff = Teff*u.K
        self.lambda_min = lambda_min*u.m
        self.lambda_max = lambda_max*u.m
        self.N = N
        
        # Calculate surface area of sun for later use
        self.sun_SA = 4*np.pi*const.R_sun**2
    
    # Function to calculate main part Planck function at given wavelength
    def Planck(self,wavelength,units=False):
        # Inputs:
        #   wavelength: wavelength to calculate Plank function at
        #   units: boolean to decide if quantities should have units
        #           (no units is preferable if using func. to integrate)
        # Returns:
        #   B: value of Planck function at that wavelength
        
        # Calculate 2hc^2 (prefactor in Planck's function)
        prefactor = (2*const.h*const.c**2)
        
        # Calculate hc/kT (constant in exponential of Planck's function)
        exp_factor = (const.h*const.c/(const.k_B*self.Teff))
        
        # Print units of B(lambda)
        #print((prefactor*wavelength**(-5)).unit)
        
        if units == False:
            # Calculate value of Planck function at this wavelength
            B = prefactor.value*wavelength.value**(-5)/(np.exp(exp_factor.value/wavelength.value)-1)
        else:
            B = prefactor*wavelength**(-5)/(np.exp(exp_factor/wavelength)-1)
            
        return(B)
    
    # Function to plot spectral energy distribution of star
    def SEDStar(self,plot=False):
        # Inputs:
        #   plot: boolean to decide to make plot of SED
        # Returns:
        #   Plot of star's SED
        #   Luminosity (integral of blackbody)*4*pi^2*R_star^2
        
        # Determine when function began running
        start_time = time.time()
        
        # Define list of wavelengths
        x = np.linspace(self.lambda_min,self.lambda_max,self.N)
        
        # Calculate Planck function at each wavelength
        y = self.Planck(x,units=True)

        # Calculate wavelength at which blackbody peaks
        lambda_peak = (2.89*(10**3)*u.micron*u.K)/self.Teff
        lambda_peak = lambda_peak.to(u.m).value
    
        # Find where blackbody peaks from my calculations
        peak_loc = np.argmax(y)
        lambda_max = x[peak_loc]
        
        # Decide whether to plot SED or not
        if plot == True:
            
            # Create wide figure (w=8in,h=4in)
            plt.figure(figsize=(8,4))
            
            # Plot data with x-axis on log scale
            plt.plot(x,y)
            plt.xscale('log')
            
            # Plot analytical and numerical wavelength peaks
            plt.vlines(lambda_peak,min(y).value,max(y).value,colors='red',\
                       linestyles='dashed',label=r"(Analytical) $\lambda_{peak}=%.2e m$"%lambda_peak)
            plt.vlines(lambda_max.value,min(y).value,max(y).value,colors='red',\
                       linestyles='dashed',label=r'(Numerical) $\lambda_{max}=%.2e m$'%lambda_max.value)
            plt.legend()
            
            # Axes labels and titles
            plt.xlabel('Wavelength ({0:latex_inline})'.format(x.unit),fontsize=14)
            plt.ylabel('Spectral Radiance ({0:latex_inline})'.format(y.unit),fontsize=14)
            plt.title(r'Spectral Energy Distribution for $T_{eff}=$%dK Star'%self.Teff.value,fontsize=18)
            plt.tight_layout()
        
        # Calculate stellar luminosity (np.trapz integrates Planck func.)
        luminosity = np.pi*self.sun_SA*np.trapz(y,x)
        print("Luminosity = {0:.3f} Lsun".format(luminosity.to(u.erg/u.s)/const.L_sun.to(u.erg/u.s)))
        
        # Calculate fraction of luminosity from below peak wavelength
        lumin_before = np.pi*self.sun_SA*np.trapz(y[:peak_loc],x[:peak_loc])
        frac_before = lumin_before/luminosity*100
        print("{0:.2f}% of energy emitted below peak".format(frac_before))
        
        # Calculate fraction of luminosity from peak wavelength and beyond
        lumin_after = np.pi*self.sun_SA*np.trapz(y[peak_loc:],x[peak_loc:])
        frac_after = lumin_after/luminosity*100
        print("{:.2f}% of energy emitted above peak \n".format(frac_after))
        
        # Tell user how long function took to run
        print('Program took %.2f sec to run'%(time.time()-start_time))
        
        return(luminosity)

test = SED()
#test.Planck(10**-6*u.m,units=False)
test.SEDStar(False)