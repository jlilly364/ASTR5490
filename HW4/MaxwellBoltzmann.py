#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 15:57:55 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import time
from scipy.integrate import quad

# Class to compute and integrate Maxwell-Boltzmann distributions
class MaxwellBoltzmann:
    
    def __init__(self,mass,temp,chemical):
        # Inputs:
        #   mass: mass of atom (in amu)
        #   temp: temperature of gas (in K)
        #   chemcial: name of atom of interest (i.e.'He','N')
        
        # Define global variables
        self.mass = mass*const.u
        self.temp = temp*u.K
        self.chemical = chemical
        
        # Calculate Earth escape speed
        self.v_esc_Earth = np.sqrt(2*const.G*const.M_earth/const.R_earth)

        # Calculate Jupiter escape speed
        self.v_esc_Jupiter = np.sqrt(2*const.G*const.M_jup/const.R_jup)
    
    # Function to calculate probability density of diff. particle speeds
    def MB(self,speed,units=False):
        # Inputs:
        #   speed: speed of particle
        # Returns:
        #   probability density of finding particle at this speed
        
        speed = speed*u.m/u.s
        
        # Calculate fraction term
        fraction = (self.mass/(2*np.pi*const.k_B.decompose()*self.temp))**(3/2)

        # Calculate exponent
        exponent = (-self.mass*(speed**2))/(2*const.k_B.decompose()*self.temp)

        # Calculate probability density at this velocity
        prob_dens = fraction*(4*np.pi*speed**2)*np.exp(exponent)
        
        if units == False:
            prob_dens = prob_dens.value
        
        return(prob_dens)
    
    # Function to plot MB distribution
    def MBDistribution(self,planet):
        # Returns:
        #   List of speeds and probability density at each speed
        
        # Determine when function began running
        start_time = time.time()
        
        # Create list of speeds
        speeds = np.linspace(10**-3,10**4,10**5)
        
        # Calculate probability density at each speed
        prob_densities = self.MB(speeds,units=True)

        # Plot MB distribution
        plt.plot(speeds,prob_densities,label='T={0}'.format(self.temp))
        plt.xlabel(r'$v\:(m\,s^{-1})$')
        plt.ylabel('Probability Density ({0:latex_inline})'.format(prob_densities.unit))
        plt.title('Maxwell-Boltzmann Distribution for {0}'.format(self.chemical))
        plt.tight_layout()
        plt.legend()
        
        # Tell user how long function took to run
        print('Program took %.2f sec to run'%(time.time()-start_time))
        
        # Choose which v_esc to use
        if planet == 'Earth':
            v_esc = self.v_esc_Earth.value
        elif planet == 'Jupiter':
            v_esc = self.v_esc_Jupiter.value
            
        # Calculate fraction of atoms with speed > v_esc for Earth
        integral,error = quad(self.MB,v_esc,np.inf)
        percentage = integral/1.0*100
        #print('{0:.2e}% of {1} atoms exceed escape velocity at T={2}'.format(integral,self.chemical,self.temp))
        
        return(speeds,prob_densities)
"""
test = MaxwellBoltzmann(4.0,2000.0,'He')
#test.MB(1000)
test.MBDistribution('Jupiter')
"""
        