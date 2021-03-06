#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:15:04 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt

# Class to generate intensity-weighted stellar profile
class LimbDarkening():
    
    def __init__(self,a,b):
        # Inputs:
        #   a: first parameter in quadratic limb darkening equation
        #   b: second parameter in quadratic limb darkening equation
        self.a = a
        self.b = b
        
    # Function to calculate quadratic limb darkening profile
    def QuadIntensity(self,x,y):
        # Inputs:
        #   x: x-coordinate to calculate intensity at
        #   y: y-coordinate to calculate intensity at
        # Returns:
        #   Intensity at that location
        
        # Set intensity at center of star
        R_star = 1.0
        
        # Calculate distance from center
        r = np.sqrt(x**2+y**2)
    
        # Calculate mu and terms that use mu
        mu = np.sqrt(1-abs(r**2/R_star**2))
    
        first_term = self.a*(1.0-mu)
        second_term = self.b*(1.0-mu)**2
        
        # Calculate intensity at r
        intensity = 1.0*(1-first_term-second_term)
       
        return(r,intensity)
    
    # Function to plot intensity of points on stellar disc
    def StarVisualize(self,gridsize,plot=True):
        # Inputs:
        #   gridsize: square dimension of grid (gridsize x gridsize)
        #   plot: boolean to choose to plot star or not
        # Returns:
        #   intensity colormap of star
        #   grid of x and y coordinates & intensities at each coordinate
        
        # Set lower bounds and size of grid
        x0, y0 = -1.01,-1.0

        # Generate list of x and y coordinates from near center to 1 R*
        x_list = np.linspace(x0,1.0,gridsize)
        y_list = np.linspace(y0,1.0,gridsize)
        x,y = np.meshgrid(x_list,y_list)
        
        # Calculate intensity at each x,y pair
        radii, intensities = self.QuadIntensity(x,y)

        print(intensities)
        if plot == True:
            # Plot color grid of intensities at each location
            plt.pcolor(x,y,intensities,shading='auto')
            cbar = plt.colorbar()
            cbar.set_label('Surface Brightness')
            plt.xlabel(r'x ($R_{star}$)')
            plt.ylabel(r'y ($R_{star}$)')
            plt.title(r'Surface Brightness of G2V at $5000\AA$')
            plt.xlim(-1.2,1.2)
            plt.ylim(-1.1,1.1)
            
        return(x,y,intensities)

star = LimbDarkening(633.27/1000,159.56/1000)
x,y,intensities = star.StarVisualize(1000,plot=False)