#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:32:14 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt

# Function to overplot transit light curves for different b's and star temps
def LightCurveCompare(rad_planet=0.05,b_values=[0.0,0.5,0.9],temperatures=['3600','5500','10000']):
    # Inputs:
    #   b_values: impact parameters of data you simulated
    #   temperatures: temperatures of stars you simulated
    # Establish figure and axes objects
    fig, (ax1,ax2,ax3) = plt.subplots(ncols=len(b_values), figsize=(15,5))
    
    # Set up empty lists of x and y data to plot
    xdata = []
    ydata = []

    # Loop over all impact parameter and star temp combinations
    for i in range(len(b_values)):
        for j in range(len(b_values)):
            # Define filename to extract data from
            filename = 'C:/Users/Jimmy/ASTR5490/HW3/TransitData/Transit_0.05Rstar_b={0}_{1}K.dat'.format(str(b_values[i]),temperatures[j])
            
            # Extract data from file and assign to x and intens
            data = np.loadtxt(filename,skiprows=1,unpack=True)
            x, intens = data[-2],data[-1]

            # Add data lists to empty list (makes a list of lists)
            xdata.append(x)
            ydata.append(intens)
    
    # Very inelegantly plot all data combinations in row of subplots
    # Loop over b=0 data
    for i in range(3):
        ax1.scatter(xdata[i],ydata[i],label='T={0}K'.format(temperatures[i]))
        j += 1
    ax1.set_ylim(0.9965,1)
    ax1.legend()
    ax1.set_xlabel(r'Horizontal Distance from Star Center ($R_{star}$)')
    ax1.set_ylabel('Relative Intensity')
    ax1.set_title('Transit of {0}'.format(rad_planet)+r'$R_{star}$ Planet'\
              +'\n'r'(b = {0})'.format(b_values[0]))
    # Loop over b=0 data
    for i in range(3,6):
        ax2.scatter(xdata[i],ydata[i],label='T={0}K'.format(temperatures[i-3]))
    ax2.set_ylim(0.9965,1)
    ax2.legend()
    ax2.set_xlabel(r'Horizontal Distance from Star Center ($R_{star}$)')
    ax2.set_ylabel('Relative Intensity')
    ax2.set_title('Transit of {0}'.format(rad_planet)+r'$R_{star}$ Planet'\
              +'\n'r'(b = {0})'.format(b_values[1]))
    # Loop over b=0 data
    for i in range(6,9):
        ax3.scatter(xdata[i],ydata[i],label='T={0}K'.format(temperatures[i-6]))
    ax3.set_ylim(0.9965,1)
    ax3.legend()
    ax3.set_xlabel(r'Horizontal Distance from Star Center ($R_{star}$)')
    ax3.set_ylabel('Relative Intensity')
    ax3.set_title('Transit of {0}'.format(rad_planet)+r'$R_{star}$ Planet'\
              +'\n'r'(b = {0})'.format(b_values[2]))
    
    # Give the plots some room to breathe
    plt.tight_layout()