#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 18:30:04 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u, constants as const
from MathTools import Planck
import time

wavelengths = np.linspace(10**-10,10**-4,10**4)*u.m
frequencies = const.c/wavelengths

sed_disks = []

for i in range()