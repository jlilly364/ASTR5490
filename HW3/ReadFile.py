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
def Read(filename,skip=2,unpack_bool=True):
    # Use numpy's loadtxt function to read columns of text file
    col1, col2, col3 = np.loadtxt(filename,skiprows=skip,unpack=unpack_bool,usecols=[0,1,2])
    return(col1,col2,col3)