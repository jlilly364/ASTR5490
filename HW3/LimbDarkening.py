#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:15:04 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
import time

# Figure out time code started to be used
start_time = time.time()

# Class to generate intensity-weighted stellar profile
class LimbDarkening():
    
    def __init__(self,star_temp,gridsize):
        # Inputs:
        #   star_temp: surface temperature of star
        #   gridsize: length and width of grid of points
        
        self.star_temp=star_temp
        self.gridsize = gridsize
        
        #   a: first parameter in quadratic limb darkening equation
        #   b: second parameter in quadratic limb darkening equation
        if self.star_temp == 5500:
            self.a = 633.27/1000
            self.b = 159.56/1000
        elif self.star_temp == 10000:
            self.a = .2481
            self.b = .2739
        elif self.star_temp == 3600:
            self.a = .626
            self.b = .226
        
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
       
        return(intensity)
    
    # Function to plot intensity of points on stellar disc
    def Star(self,plot=True):
        # Inputs:
        #   plot: boolean to choose to plot star or not
        # Returns:
        #   intensity colormap of star
        #   grid of x and y coordinates & intensities at each coordinate
        
        # Set lower bounds and size of grid
        x0, y0 = -1.0,-1.0

        # Generate list of x and y coordinates from near center to 1 R*
        x_list = np.linspace(x0,1.0,self.gridsize)
        y_list = np.linspace(y0,1.0,self.gridsize)
        x,y = np.meshgrid(x_list,y_list)
        
        # Calculate intensity at each x,y pair
        intensities = self.QuadIntensity(x,y)

        #print(intensities)
        if plot == True:
            # Plot color grid of intensities at each location
            plt.pcolor(x,y,intensities)#,shading='nearest')
            cbar = plt.colorbar()
            cbar.set_label('Surface Brightness',fontsize=14)
            plt.xlabel(r'x ($R_{star}$)',fontsize=14)
            plt.ylabel(r'y ($R_{star}$)',fontsize=14)
            plt.title('Surface Brightness of T={0}K Star at '.format(self.star_temp)+r'$5000\AA$',fontsize=18)
            plt.xlim(-1.2,1.2)
            plt.ylim(-1.1,1.1)
            
        return(x,y,intensities)
    
    # Place star at particular point in 
    def Transit(self,rad_planet,b,plot=False):
        
        self.b = b
        
        # Generate intensity-weighted coordinate grid
        x_grid,y_grid,intensities_star = self.Star(plot=False)
        real_intensities_star = [z for z in intensities_star.flatten() if ~np.isnan(z)]
        original_total = np.sum(real_intensities_star)
        
        # Make empty list of relative intensities
        light_curve = []
        
        # Initialize figure and axis object for plotting
        fig, ax = plt.subplots()
        
        # Set loop counter
        i=0
        
        # Establish array to add data to
        data = np.zeros([len(x_grid[0]), 5])
        
        # Calculate intensity from visible star throughout transit
        for x in x_grid[0]:
            
            # Identify location of planet center
            planet_center = [-x,self.b]
            
            # Calculate x,y, and total distances of all points from planet center
            xdist = x_grid - planet_center[0]
            ydist = y_grid - planet_center[1]
            dist = np.sqrt(xdist**2+ydist**2)
            
            # Find pixels within planet radius
            planet_ids = np.where(dist < rad_planet)
            
            # Set intensity to 0 wherever planet blocks the star        
            non_transit_intensities = intensities_star[planet_ids]
            intensities_star[planet_ids] = 0
            
            # Remove NaNs from intensity list
            real_intensities_transit = [z for z in intensities_star.flatten() if ~np.isnan(z)]
            
            # Calculate total observed intensity at this point in transit
            transit_total = np.sum(real_intensities_transit)
            
            # Calculate relative intensity to non-transit
            light_fraction = transit_total/original_total
            light_curve.append(light_fraction)
            
            # Print status of loop
            print("Now completing Loop {0} out of {1}: Rel. Intens. = {2:.2f}".format(i,len(x_grid[0]),light_fraction))

            # Plot star with planet in front if user desires
            if plot==True:
                fig.clf()
                ax.pcolor(x_grid,y_grid,intensities_star,cmap='hot')#,shading='nearest')
                ax.set_xlim(-1.2,1.2)
                ax.set_ylim(-1.1,1.1)
                ax.set_facecolor('black')
                #fig.savefig("C:/Users/Jimmy/Downloads/Test/test_{0}.png".format(i),)
                fig.canvas.draw()
                plt.close(fig)
            
            # Save important data
            data[i] = self.star_temp, rad_planet, self.b, -x, light_fraction
            
            # Reset intensities to original
            intensities_star[planet_ids] = non_transit_intensities    
            i += 1
            
            # Determine how long program has been running
            """looptime = time.time() - start_time
            print("Time elapsed: {0:.3f}".format(looptime))"""
        
        # Save important data to text file (only has to be run once)
        fileout = 'C:/Users/Jimmy/ASTR5490/HW3/TransitData/Transit_{0}Rstar_b={1}_{2}K.dat'.format(rad_planet,self.b,self.star_temp)
        np.savetxt(fileout, data, fmt = "%11.2f %11.2f %11.2f %11.9f %11.9f",comments='#',
               header="{:^10s}{:^11s}{:^11s}{:^11s}{:^11s}"\
                      .format('star_temp(K)','rad_planet(R*)', 'b', 'x_pos', 'rel_intens'))
        
        # Plot transit light curve
        plt.scatter(x_grid[0],light_curve)
        plt.xlabel(r'Horizontal Distance from Star Center ($R_{star}$)',fontsize=14)
        plt.ylabel('Relative Intensity',fontsize=14)
        plt.title('Transit of {0}'.format(rad_planet)+r'$R_{star}$ Planet'\
                  +'\n'+r'($T_{star}$ = '+'{0}K, b = {1})'.format(self.star_temp,self.b),fontsize=18)
        #plt.savefig('Transit_{0}Rstar_b={1}_{2}K.png'.format(rad_planet,self.b,self.star_temp))
        
        # Determine how long it took the program to run
        runtime = time.time() - start_time
        print("My program took {0:.2f} seconds to run".format(runtime))