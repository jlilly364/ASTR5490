#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:29:38 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt

# Class to simulate tranit of planet in front of star
class Transit():
    
    def __init__(self,rad_planet,b):
        # Inputs:
        #   rad_planet: radius of transiting planet (fraction of R_star)
        #   b: impact parameter of planet relative to star
        
        # Define global variables
        self.Rp = rad_planet
        self.b = b
        
    # Function to calculate quadratic limb darkening profile
    def PlanetIntensity(self,x,y):
        # Inputs:
        #   x: x-coordinate to calculate intensity at
        #   y: y-coordinate to calculate intensity at
        # Returns:
        #   Intensity at that location
        
        # Calculate distance from center
        r = np.sqrt(x**2+y**2)

        # Calculate intensity at r
        intensity = np.ones((len(x),len(y)))

        return(intensity)
        
    def Planet(self,gridsize,buffer=1.25):
        # Set bounds and size of coordinate grid
        xmin, ymin = -self.Rp*buffer,-self.Rp*buffer
        xmax, ymax = self.Rp*buffer,self.Rp*buffer

        # Generate list of x and y coordinates from near center to 1 R*
        x_list = np.linspace(xmin,xmax,gridsize)
        y_list = np.linspace(ymin,ymax,gridsize)
        x,y = np.meshgrid(x_list,y_list)
        
        distance = np.sqrt(x**2+y**2)
        
        inside = distance <= self.Rp
        
        planet_intensities = self.PlanetIntensity(x[inside], y[inside])
        #print(len(x),len(y),len(planet_intensities))
        
        #inside = planet_radii < self.Rp
        #print(len(x),len(y),len(planet_intensities))
        #print(planet_intensities.shape)
        
        #inside = np.where(planet_radii > self.Rp)
        
        #print(len(x[inside]),len(y[inside]),len(intensities[inside]))
        
        #plt.pcolor(x,y,planet_intensities,shading='nearest')
        plt.pcolor(x[inside],y[inside],planet_intensities,shading='nearest')
        plt.xlim(-.25,.25)
        plt.ylim(-.25,.25)
        
test = Transit(0.05,0)
test.Planet(100)