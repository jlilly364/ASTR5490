#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:32:14 2020

@author: jimmy
"""
x_values = []
intensities = []
temps = ['3600','5500','10000']

def LightCurveCompare(b,temperatures):
    x_values = []
    intensities = []
    for temp in temperatures:
        filename = 'Transit_0.05Rstar_b={0}_{1}K.dat'.format(str(b),temp)
        data = np.loadtxt(filename,skiprows=1,unpack=True)
        x, intens = data[-2],data[-1]
        x_values.append(x)
        intensities.append(intens)
    #plt.scatter(x_values[0],intensities[0])#,label=temp+'K')
    print(intensities[0])

LightCurveCompare(0.0,temps)