#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:16:09 2020

@author: jimmy
"""
# Import numpy module
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

# Function to read simple text files
def Read(filename,col_names,unpack_bool=True):
    # Use numpy's loadtxt function to read columns of text file
    data = np.genfromtxt(filename,delimiter=' ',names=True)
    col1, col2, col3 = data[col_names[0]],data[col_names[1]],data[col_names[2]]
    return(col1,col2,col3)

# Function to read NASA Exoplanet Archive text files
def ReadNASA(filename,skip):
    # Read and return data of confirmed exoplanets from NASA Exoplanet Archive
    data = np.genfromtxt(filename,dtype=None,delimiter=',',skip_header=skip,names=True,invalid_raise=False,encoding=None)
    return data

# Function to read HW4 bandpass data (Kepler and Spitzer 4.5um)
def ReadBandpass(filename,xvariable,delim='\t',skip=8):
    
    # Extract data from file
    data = np.genfromtxt(filename,delimiter=delim,skip_header=skip)
    
    # Extract data into 2 separate lists
    xdata, ydata = zip(*data)
    
    # Save lists as numpy arrays (convert first list to m from um)
    xdata, ydata = np.asarray(xdata)/10**6*u.m, np.asarray(ydata) # wavelengths in micron
    
    # Convert wavelengths to frequencies
    if xvariable == 'freq':
        xdata = xdata.to(u.s**(-1), equivalencies=u.spectral())
    
    return(xdata,ydata)