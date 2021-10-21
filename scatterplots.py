#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 09:23:45 2021

@author: fra
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results = pd.read_csv( 'phylogenetic_community1.txt',sep=',')


stats = pd.read_csv('PARAMS_COMMUNITY1_ACCEPTED/Outputs/statistics_runs.csv')

stats.replace([np.inf, -np.inf], np.nan, inplace=True)
stats.fillna(0,inplace=True)



for k in [1,5,9,13,17]:
    fig, axs = plt.subplots(2,2)
    plt.subplots_adjust(left=0.12,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.4)
    axs = axs.flatten()
  
    for i in range(k,k+4):
      j=i%4
      print(j)
      #minx =min( min(stats.values[:,i]),results.values[0][i])
      #maxx = max (max(stats.values[:,i]),results.values[0,i])
      
      axs[j].hist(stats.values[:,i])
      axs[j].set_ylabel(stats.keys()[i])
      axs[j].axvline(x=results.values[0,i], linestyle="dashed", color="red")
      #axs[j].set_ylim( (minx,maxx) )