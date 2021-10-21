#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 16:43:31 2021

@author: fra
"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

params_sex=pd.read_csv('PARAMS_vary_cmulti/param_processed_patch0_partnerships.csv',sep=' ')

numberevents=12

assortativity=params_sex.assortativity[:numberevents]
c_multi=params_sex.c_multiplier[:numberevents]
breakup_scale=params_sex.breakup_scale_multiplier_overall[:numberevents]



mean=np.zeros(numberevents)
var=np.zeros(numberevents)
median=np.zeros(numberevents)
twofive=np.zeros(numberevents)
twentyfive=np.zeros(numberevents)
seventyfive=np.zeros(numberevents)
ninentyseven=np.zeros(numberevents)
var_offspring_distribution=np.zeros(numberevents)

plt.figure()    

for i in range(1,numberevents+1):

    stats = pd.read_csv('PARAMS_vary_cmulti/Outputs/output_Run%d_transmission_network.tsv'%i,sep='\t', names=['Infector','Infected','Time'])
    
    condition1 = stats.Time>1980
    condition2 = stats.Time<2011 #Right data censoring 
    
    stats1=stats[condition1]
    #Stop the count to those infected prior to 2011
    potentialinfectors=stats[condition2].iloc[:,0].values

    stats=stats1[stats.Infector.isin(potentialinfectors)]

    
    grouped=stats1.groupby(by='Infector').size().to_list()
    
    numberofzeros=len(stats[condition1].Infector) - len(grouped) #(everyone who is not included did not infected anybody)
    
    for j in range(numberofzeros):
        grouped.append(0)
    '''
    print(max(grouped))
    if i==1 or i==5 or i==10:
        plt.hist(grouped,label = 'assortativity = %.2f'%c_multi[i])
    plt.legend()
    
    plt.title('assortativity: %.2f'%assortativity[i-1])
    '''
    twofive[i-1],twentyfive[i-1],median[i-1],seventyfive[i-1],ninentyseven[i-1]=np.quantile(grouped,[0.025,0.25,0.5,0.75,0.975])
    mean[i-1] = np.mean(grouped)
    var[i-1] = np.var(grouped)
    

fig, (ax1, ax2) = plt.subplots(2)

m1,b1=np.polyfit(c_multi,mean,1)
m2,b2=np.polyfit(c_multi,var,1)

c_m=np.sort(c_multi)
ax1.scatter(c_m,mean, color='orange', label='offspring distribution average size')
ax2.scatter(c_m, var, color='blue', label='offspring distribution variance')
#ax1.plot(c_m,m1*c_m+b1, label='m=%.2f'%m1)
#ax2.plot(c_m,m2*c_m+b2, label='m=%.2f'%m2)


ax1.set_xlabel('missreporting')
ax2.set_xlabel('missreporting')
ax1.set_ylabel('N')
ax1.legend()
ax2.legend()


fig, (ax1, ax2) = plt.subplots(2)


m1,b1=np.polyfit(assortativity,mean,1)
m2,b2=np.polyfit(assortativity,var,1)

ass=np.sort(assortativity)
ass=np.sort(assortativity)
ax1.scatter(assortativity,mean, color='orange', label='offspring distribution average size')
ax2.scatter(assortativity, var, color='blue', label='offspring distribution variance')
ax1.plot(ass,m1*ass+b1, label='m=%.2f'%m1)
ax2.plot(ass,m2*ass+b2, label='m=%.2f'%m2)

ax1.set_xlabel('assortativity')
ax2.set_xlabel('assortativity')
ax1.set_ylabel('N')
ax1.legend()
ax2.legend()

'''
plt.figure()
plt.scatter(c_multi,assortativity)



stats_community1 = pd.read_csv('PARAMS_vary_assortativity/Outputs/statistics_runs.csv',sep=',')

width=0.35
meanstats = stats.mean().to_list()
fig,ax = plt.subplots()
x= np.arange(len(meanstats)-2)
ax.scatter(x-width/2, meanstats[:-2],color='red')
ax.scatter(x-width/2, stats_community1,color='blue')
labels=['run','max_H','min_H','a_BL_mean','a_BL_median','a_BL_var','i_BL_mean_i','i_BL_var_i','i_BL_mean_e','i_BL_var_e','colless','sackin','WD_ratio','DeltaW',
           'max_ladder','staircaseness_1','staircaseness_2','max_L','t_max_L','slope_1']
stats= pd.read_csv('PARAMS_vary_assortativity/Outputs/statistics_runs.csv',sep=',', names=labels)


fig, axs = plt.subplots(4)
for index,k in enumerate(stats.columns[16:21]):
    axs[index].scatter(assortativity,stats[k])
    axs[index].set_title(stats.keys()[index+16])
    if index==4:
        axs[index].set_xlabel('assortativity')




ax.set_xticks(x-width/2)
plt.xticks(rotation='vertical')
ax.set_xticklabels(labels)
#Topology
fig, axs = plt.subplots(5)
for index,k in enumerate(stats.columns[1:6]):
    axs[index].scatter(assortativity,stats[k])
    axs[index].set_title(stats.keys()[index+1])
    if index==4:
        axs[index].set_xlabel('assortativity')

#Branch length
fig, axs = plt.subplots(4)
for index,k in enumerate(stats.columns[:4]):
    axs[index].scatter(assortativity,stats[k])
    axs[index].set_title(stats.keys()[index+4])
    if index==4:
        axs[index].set_xlabel('assortativity')
        
'''