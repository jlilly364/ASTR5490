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
    
    def __init__(self,xvariable,yvariable,Teff=5780,lambda_min=10**-10,lambda_max=10**-4,N=10**5):
        # Inputs:
        #   xvariable: 'freq' or 'wavelen' to determine which Planck form to use
        #   yvariable: 'planck' or 'luminosity' or 'xvar_lumin'
        #   T_eff: effective temperature of blackbody (in K)
        #   lambda_min: minimum wavelength to calc. Planck function over (in m)
        #   lambda_max: minimum wavelength to calc. Planck function over (in m)
        #   Number of subintervals to integrate over    
    
        # Cast initial parameters as global variables
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.Teff = Teff*u.K
        self.lambda_min = lambda_min*u.m
        self.lambda_max = lambda_max*u.m
        self.N = N
        
        # Calculate wavelength at which blackbody peaks (in m)
        self.lambda_peak = ((2.90*(10**3)*u.micron*u.K)/self.Teff).to(u.m)
        
        # Set boundaries of analysis depending on x variable
        if self.xvariable == 'wavelen':
            self.x_min = self.lambda_min
            self.x_max = self.lambda_max
            
            # Calculate wavelength at which blackbody peaks (in m)
            self.y_peak = self.lambda_peak
            
            # Define x-axis label
            self.xlabel = r'$\lambda$ ({0:latex_inline})'.format(self.x_min.unit)
            
            # Define v-line labels
            self.vlabel = r'Analytical (Planck): $\lambda_{{max}}={0:.2e}$ {1:latex_inline}'.format(self.y_peak.value,self.y_peak.unit)
            self.num_vlabel = r'Numerical: $\lambda_{{max}}={0:.2e}$ {1:latex_inline}'
            
            if self.yvariable == 'planck':
                self.ylabel = r'Spectral Radiance per $\Omega$ ({0:latex_inline})'
            elif self.yvariable == 'luminosity':
                self.ylabel = r'$L_{{\lambda}}$ ({0:latex_inline})'
            elif self.yvariable == 'xvar_lumin':
                self.ylabel = r'$\lambda L_{{\lambda}}$ ({0:latex_inline})'
        
        elif self.xvariable == 'freq':
            # Convert wavelength range to frequency range
            self.x_min = self.lambda_min.to(u.s**(-1), equivalencies=u.spectral())
            self.x_max = self.lambda_max.to(u.s**(-1), equivalencies=u.spectral())
            
            # Calculate frequency at which blackbody peaks (in s^-1)
            self.y_peak = (5.88*(10**10)*(u.s**(-1))/u.K)*self.Teff
            
            # Define x-axis label
            self.xlabel = r'$\nu$ ({0:latex_inline})'.format(self.x_min.unit)
            
            # Define v-line labels
            self.vlabel = r'Analytical (Planck): $\nu_{{max}}={0:.2e}$ {1:latex_inline}'.format(self.y_peak.value,self.y_peak.unit)
            self.num_vlabel = r'Numerical: $\nu_{{max}}={0:.2e}$ {1:latex_inline}'
            
            if self.yvariable == 'planck':
                self.ylabel = r'Spectral Radiance per $\Omega$ ({0:latex_inline})'
            elif self.yvariable == 'luminosity':
                self.ylabel = r'$L_{{\nu}}$ ({0:latex_inline})'
            elif self.yvariable == 'xvar_luminos':
                self.ylabel = r'$\nu L_{{\nu}}$ ({0:latex_inline})'
        else:
            print("Valid entries are 'wavelen' or 'freq'")
        
        # Calculate surface area of sun for later use
        self.sun_SA = 4*np.pi*const.R_sun**2
    
    # Function to calculate main part Planck function at given wavelength
    def Planck(self,x,units=True):
        # Inputs:
        #   x: value of x-variable to calculate Planck function at
        #   units: boolean to decide if quantities should have units
        #           (no units is preferable if using func. to integrate)
        # Returns:
        #   B: value of Planck function at that wavelength
        
        if self.xvariable == 'wavelen':
            # Calculate 2hc^2 (prefactor in Planck's function)
            prefactor = (2*const.h*const.c**2)
            
            # Calculate hc/kT (constant in exponential of Planck's function)
            exp_factor = (const.h*const.c/(const.k_B*self.Teff))
            
            if units == False:
                # Calculate value of Planck function at this wavelength
                B = prefactor.value*x.value**(-5)/(np.exp(exp_factor.value/x.value)-1)
            else:
                B = prefactor*x**(-5)/(np.exp(exp_factor/x)-1)
        
        elif self.xvariable == 'freq':
            
            # Calculate 2h/c^2 (prefactor in Planck's function)
            prefactor = 2*const.h/const.c**2
            
            # Calculate h/kT (constant in exponential of Planck's function)
            exp_factor = const.h/(const.k_B*self.Teff)
        
            if units == False:
                # Calculate value of Planck function at this wavelength
                B = prefactor.value*x.value**3/(np.exp(exp_factor.value*x.value)-1)
            else:
                B = prefactor*x**3/(np.exp(exp_factor*x)-1)
            
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
        if self.xvariable == 'wavelen':
            x = np.linspace(self.x_min,self.x_max,self.N)
        elif self.xvariable == 'freq':
            x = np.linspace(self.x_max,self.x_min,self.N)
            
        # Calculate Planck function at each wavelength
        y = self.Planck(x,units=True)
        
        if self.yvariable == 'luminosity':
            # Convert Planck function to luminosity
            y *= np.pi*self.sun_SA
            y = np.multiply(x,y)
        elif self.yvariable == 'xvar_luminos':
            # Convert Planck function to luminosity*xvariable (planck * x**2)
            y *= np.pi*self.sun_SA
            y = np.multiply(np.square(x),y)
            
        # Find where blackbody peaks from my calculations
        peak_loc = np.argmax(y)
        numerical_max = x[peak_loc]
        
        if self.yvariable == 'planck':
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
        
        # Decide whether to plot SED or not
        if plot == True:
            
            # Create wide figure (w=8in,h=4in)
            plt.figure(figsize=(8,4))
            
            # Plot data with x-axis on log scale
            plt.plot(x,y)
            plt.xscale('log')
            
            # Plot analytical and numerical wavelength peaks
            plt.vlines(self.y_peak.value,min(y).value,max(y).value,colors='green',\
                       linestyles='dashed',label=self.vlabel)
            plt.vlines(numerical_max.value,min(y).value,max(y).value,colors='red',\
                       linestyles='dashed',label=self.num_vlabel.format(numerical_max.value,numerical_max.unit))
            plt.legend()
            
            # Axes labels and titles
            plt.xlabel(self.xlabel,fontsize=14)
            plt.ylabel(self.ylabel.format(y.unit),fontsize=14)
            plt.title(r'Spectral Energy Distribution for $T_{eff}=$%dK Star'%self.Teff.value,fontsize=18)
            plt.tight_layout()
        
        # Tell user how long function took to run
        print('Program took %.2f sec to run'%(time.time()-start_time))
        
        return(x,y)

# Code to test class and functions
"""test = SED('freq','xvar_luminos')
#test.Planck(10**-6*u.m,units=False)
test.SEDStar(True)"""