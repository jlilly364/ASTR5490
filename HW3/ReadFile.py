#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:16:09 2020

@author: jimmy
"""
# Import numpy module
import numpy as np
import matplotlib.pyplot as plt

# Function to read simple text files
def Read(filename,col_names,unpack_bool=True):
    # Use numpy's loadtxt function to read columns of text file
    data = np.genfromtxt(filename,delimiter='\t',names=True)
    col1, col2, col3 = data[col_names[0]],data[col_names[1]],data[col_names[2]]
    return(col1,col2,col3)

# Function to read NASA Exoplanet Archive text files
def ReadNASA(filename,skip):
    # Read and return data of confirmed exoplanets from NASA Exoplanet Archive
    data = np.genfromtxt(filename,dtype=None,delimiter=',',skip_header=skip,names=True,invalid_raise=False,encoding=None)
    return data