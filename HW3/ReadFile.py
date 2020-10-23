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
def Read(filename,unpack_bool=True):
    # Use numpy's loadtxt function to read columns of text file
    data = np.genfromtxt(filename,delimiter=' ',names=True)
    col1, col2, col3 = data['col0'],data['col1'],data['col2']
    return(col1,col2,col3)