##############################################################################
# Figure 1
##############################################################################
import numpy as np
import pandas as pd
import scipy as sp
from scipy import signal
from scipy import stats
from scipy.io import loadmat
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math
from sklearn import metrics
from sklearn import linear_model
import pickle
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import cbook
from matplotlib import colors
from matplotlib.colors import Normalize

# FUNCTIONS: 
    # - Fraction of trials 
    # - Opto_plots

##############################################################################
## Fraction of trials
##############################################################################
def Fraction_of_trials(mice, master_log):
    trialtypes1=['SomHit','SomCR','SomMiss','SomFA']
    trialtypes2=['VisHit',  'VisCR','VisMiss','VisFA']
    total_number_sessions=0
    temp_units1=[]
    temp_units2=[]
    trial_counter=0
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])] #you can just use day log because having multiple units will not change the fraction           
            temp=[]
            for i in np.arange(len(trialtypes1)):
                temp.append(len(day_log[day_log['Stim/Block/Response']==trialtypes1[i]].index))
            temp=[x/sum(temp) for x in temp]
            temp_units1.append(temp)
            temp=[]
            for i in np.arange(len(trialtypes2)):
                temp.append(len(day_log[day_log['Stim/Block/Response']==trialtypes2[i]].index))
            temp=[x/sum(temp) for x in temp]
            temp_units2.append(temp)
            total_number_sessions+=1
            
    
    Mean1=np.mean(temp_units1,axis=0)
    Mean2=np.mean(temp_units2,axis=0)
    SEM1=np.std(temp_units1, axis=0)/np.sqrt(len(temp_units1))
    SEM2=np.std(temp_units2, axis=0)/np.sqrt(len(temp_units2))
    
    HITs=[Mean1[0], Mean2[0]]
    HITsSEM=[SEM1[0], SEM2[0]]
    CRs=[Mean1[1], Mean2[1]]
    CRsSEM=[SEM1[1], SEM2[1]]
    MISSs=[Mean1[2], Mean2[2]]
    MISSsSEM=[SEM1[2], SEM2[2]]
    FAs=[Mean1[3], Mean2[3]]
    FAsSEM=[SEM1[3], SEM2[3]]
    
    ind = np.arange(2)    # the x locations for the groups
    width = 0.5       # the width of the bars: can also be len(x) sequence
    
    fig, ax=plt.subplots(1,1, figsize=(2,6))
    plt.sca(ax)
    p1 = plt.bar(ind, HITs, width, yerr=HITsSEM, color='k', alpha=0.9)
    p2 = plt.bar(ind, CRs, width, bottom=HITs, yerr=CRsSEM,  color='k',alpha=0.7)
    p3 = plt.bar(ind, MISSs, width, bottom= [sum(x) for x in zip(CRs, HITs)], yerr=MISSsSEM, color='k', alpha=0.4)
    p4 = plt.bar(ind, FAs, width, bottom=[sum(x) for x in zip(MISSs,CRs, HITs)], yerr=FAsSEM, color='k', alpha=0.1)
    
    plt.ylabel('Fraction of trials')
    plt.xticks(ind, ('Touch\nBlock', 'Visual\nBlock'))
    plt.xlim(-0.5,1.5)
    plt.title('n=9 mice, '+str(total_number_sessions)+' sessions')
    #plt.legend((p1[0], p2[0], p3[0],p4[0]), ('Hit', 'CR', 'Miss', 'FA'), loc='upperright')
    
    return HITs, HITsSEM, CRs, CRsSEM, MISSs, MISSsSEM, FAs, FAsSEM

##############################################################################
## Opto_plots
##############################################################################

def Opto_plots(mice, master_log, All_units_OptoMetrics_df):
    OptoTag_list=[]
    OptoTagSameTT_list=[]
    SALT_list=[]
    SALTSameTT_list=[]
    InBetween_list=[]
    Rest=[]
    count=0
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                if unit_log['Category'].values[0] == 'OptoTag':
                    OptoTag_list.append(unit_log['unit_name'].values[0])
                    count+=1
                if unit_log['Category2'].values[0] == 'OptoNetwork':
                    SALT_list.append(unit_log['unit_name'].values[0])
                    count+=1
                elif unit_log['Category'].values[0] == 'SameTT':
                    OptoTagSameTT_list.append(unit_log['unit_name'].values[0])
                    count+=1
                elif unit_log['Category3'].values[0] == 'OptoNetwork_SameTT':
                    SALTSameTT_list.append(unit_log['unit_name'].values[0])      
                    count+=1
                elif unit_log['Category4'].values[0] == 'InBetween':
                    InBetween_list.append(unit_log['unit_name'].values[0])
                    count+=1
    
                else:
                    Rest.append(unit_log['unit_name'].values[0])
    
    count=0               
    for each in All_units_OptoMetrics_df.index.values:
        if All_units_OptoMetrics_df.loc[each, 'unit_name'] in Rest:
            All_units_OptoMetrics_df.at[each,'Category']= 'Rest'
        elif All_units_OptoMetrics_df.loc[each, 'unit_name'] in OptoTag_list:
            All_units_OptoMetrics_df.at[each,'Category']= 'OptoTag'  
            count+=1
        elif All_units_OptoMetrics_df.loc[each, 'unit_name'] in SALT_list:
            All_units_OptoMetrics_df.at[each,'Category']= 'SALT'
            count+=1
        elif All_units_OptoMetrics_df.loc[each, 'unit_name'] in OptoTagSameTT_list:
            All_units_OptoMetrics_df.at[each,'Category']= 'OptoTagSameTT_list'
            count+=1
        elif All_units_OptoMetrics_df.loc[each, 'unit_name'] in SALTSameTT_list:
            All_units_OptoMetrics_df.at[each,'Category']= 'SALTSameTT_list'
            count+=1
        elif All_units_OptoMetrics_df.loc[each, 'unit_name'] in InBetween_list:
            All_units_OptoMetrics_df.at[each,'Category']= 'InBetween'
            count+=1
    
    
        else:
            All_units_OptoMetrics_df.at[each,'Category']= 'NotClaustrum'
      
    a=0
    b=0
    c=0
    d=0
    e=0
    f=0
    
    for unit in  np.unique(All_units_OptoMetrics_df['unit_name']):
        unit_log=All_units_OptoMetrics_df[All_units_OptoMetrics_df['unit_name']==unit]
        if unit_log['Category'].values[0]=='OptoTag':
            a+=1
        elif unit_log['Category'].values[0]=='SALT':
            b+=1
        elif unit_log['Category'].values[0]=='OptoTagSameTT_list':
            c+=1
        elif unit_log['Category'].values[0]=='SALTSameTT_list':
            d+=1
        elif unit_log['Category'].values[0]=='InBetween':
            e+=1
        elif unit_log['Category'].values[0]=='NotClaustrum':
            f+=1
    
    # Example not modulated (raster, mean waveform with few individual)
    Sub_df=All_units_OptoMetrics_df[All_units_OptoMetrics_df['frequency']==2]
    Sub_df=Sub_df[Sub_df['Category']=='InBetween']
    example=7
    example_df=Sub_df.loc[Sub_df.index[example], :]
    fig, ax= plt.subplots(1,1,figsize=(8,8))
    ax.fill([0,0.003, 0.003, 0],[0,0,150,150], color='dodgerblue', alpha=0.2, edgecolor='w' )
    for i,pulse in enumerate(example_df.pulses[:150]):
        pulse_aligned_spikes=example_df.all_spikes-pulse
        ax.scatter(pulse_aligned_spikes, np.zeros_like(pulse_aligned_spikes)+i, c='grey')
    ax.set_xlim(-0.02,0.05)
    ax.set_ylim(0,150)
    ax.set_xlabel('Time (ms)', size=24)
    ax.set_ylabel('Pulses', size=24)
    ax.set_xticks([-0.02, -0.01 , 0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax.set_xticklabels(['-20','-10','0','10','20','30','40', '50'], size=20)
    ax.set_yticks([0,50,100,150])
    ax.set_yticklabels(['0','50','100','150'], size=20)
    ax.set_title('Example no response')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
    
    # Example SALT (raster, mean waveform with few individual)
    Sub_df=All_units_OptoMetrics_df[All_units_OptoMetrics_df['frequency']==2]
    Sub_df=Sub_df[Sub_df['Category']=='SALT']

    example=18
    example_df=Sub_df.loc[Sub_df.index[example], :]
    fig, ax= plt.subplots(1,1,figsize=(8,8))
    ax.fill([0,0.003, 0.003, 0],[0,0,150,150], color='dodgerblue', alpha=0.2, edgecolor='w' )
    for i,pulse in enumerate(example_df.pulses[:150]):
        pulse_aligned_spikes=example_df.all_spikes-pulse
        ax.scatter(pulse_aligned_spikes, np.zeros_like(pulse_aligned_spikes)+i, c='c')
    ax.set_xlim(-0.02,0.05)
    ax.set_ylim(0,150)
    ax.set_xlabel('Time (ms)', size=24)
    ax.set_ylabel('Pulses', size=24)
    ax.set_xticks([-0.02, -0.01 , 0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax.set_xticklabels(['-20','-10','0','10','20','30','40', '50'], size=20)
    ax.set_yticks([0,50,100,150])
    ax.set_yticklabels(['0','50','100','150'], size=20)
    ax.set_title('Example SALT')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
    # Example Opto (raster, mean waveform with few individual, Summary of Corr/proba/Latency)
    Sub_df=All_units_OptoMetrics_df[All_units_OptoMetrics_df['frequency']==10]
    Sub_df=Sub_df[Sub_df['Category']=='OptoTag']
    example=7
    example_df=Sub_df.loc[Sub_df.index[example], :]
    fig, ax= plt.subplots(1,1,figsize=(8,8))
    ax.fill([0,0.003, 0.003, 0],[0,0,150,150], color='dodgerblue', alpha=0.2, edgecolor='w' )
    for i,pulse in enumerate(example_df.pulses[:150]):
        pulse_aligned_spikes=example_df.all_spikes-pulse
        ax.scatter(pulse_aligned_spikes, np.zeros_like(pulse_aligned_spikes)+i, c='mediumpurple')
    ax.set_xlim(-0.02,0.05)
    ax.set_ylim(0,150)
    ax.set_xlabel('Time (ms)', size=24)
    ax.set_ylabel('Pulses', size=24)
    ax.set_xticks([-0.02, -0.01 , 0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax.set_xticklabels(['-20','-10','0','10','20','30','40', '50'], size=20)
    ax.set_yticks([0,50,100,150])
    ax.set_yticklabels(['0','50','100','150'], size=20)
    ax.set_title('Example Optotag')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    ax=ax.flatten()
    
    plt.sca(ax[0])
    plt.plot(np.arange(len(example_df['mean_cont_waveform1'])), example_df['mean_cont_waveform1'],color='grey')
    plt.plot(np.arange(len(example_df['mean_evoked_waveform1'])), example_df['mean_evoked_waveform1'],color='dodgerblue')
    plt.ylim(-0.1,0.03)
    plt.axis('off')
    
    plt.sca(ax[1])
    plt.plot(np.arange(len(example_df['mean_cont_waveform2'])), example_df['mean_cont_waveform2'],color='grey')
    plt.plot(np.arange(len(example_df['mean_evoked_waveform2'])), example_df['mean_evoked_waveform2'],color='dodgerblue')
    plt.ylim(-0.1,0.03)
    plt.axis('off')
    
    plt.sca(ax[2])
    plt.plot(np.arange(len(example_df['mean_cont_waveform3'])), example_df['mean_cont_waveform3'],color='grey')
    plt.plot(np.arange(len(example_df['mean_evoked_waveform3'])), example_df['mean_evoked_waveform3'],color='dodgerblue')
    plt.ylim(-0.1,0.03)
    plt.xticks([0,15,30], ['0','0.5','1'])
    #plt.axis('off')
    
    plt.sca(ax[3])
    plt.plot(np.arange(len(example_df['mean_cont_waveform4'])), example_df['mean_cont_waveform4'],color='grey')
    plt.plot(np.arange(len(example_df['mean_evoked_waveform4'])), example_df['mean_evoked_waveform4'],color='dodgerblue')
    plt.axis('off')
    
    
    #plot 2d: latency vs reliability
    Sub_df=All_units_OptoMetrics_df[All_units_OptoMetrics_df['frequency']==10]
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    groups = Sub_df.groupby('Category')
    ax.set_prop_cycle(color=['grey','grey','grey','grey','grey','grey', 'grey'])
    for name, group in groups:
        ax.scatter(group.mean_latencies,group.reliability, s=20, alpha=0.5, linewidths=0)
    ax.scatter(Sub_df[Sub_df['Category']=='OptoTag'].mean_latencies, Sub_df[Sub_df['Category']=='OptoTag'].reliability, s=20, c='mediumpurple', alpha=1)
    plt.show()
    ax.set_xticks([0,0.005, 0.01, 0.015, 0.02])
    ax.set_xticklabels(['0','5','10','15', '20'], size=20)
    ax.set_yticklabels(['0','0.2','0.4','0.6', '0.8','1'], size=20)
    ax.set_xlim(0,0.020)
    ax.set_ylim(0,1)
    ax.set_xlabel('Latency (msec)', size=24)
    ax.set_ylabel('Reliability', size=24)
    ax.set_title('Summary_plot_LatencyvsReliability')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
    #plot 2d: latency vs Correlation
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    groups = Sub_df.groupby('Category')
    ax.set_prop_cycle(color=['grey','grey','grey','grey','grey','grey', 'grey'])
    for name, group in groups:
        ax.scatter(group.mean_latencies,group.waveform_corr, s=20, alpha=0.5, linewidths=0)
    ax.scatter(Sub_df[Sub_df['Category']=='OptoTag'].mean_latencies, Sub_df[Sub_df['Category']=='OptoTag'].waveform_corr, s=20, c='mediumpurple', alpha=1)
    plt.show()
    ax.set_xticks([0,0.005, 0.01, 0.015, 0.02])
    ax.set_xticklabels(['0','5','10','15', '20'], size=20)
    ax.set_yticklabels(['0','0.2','0.4','0.6', '0.8','1'], size=20)
    ax.set_xlim(0,0.020)
    ax.set_ylim(0,1.01)
    ax.set_xlabel('Latency (msec)', size=24)
    ax.set_ylabel('Waveform correlation', size=24)
    ax.set_title('Summary_plot_LatencyvsWVcorr')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return a,b,c,d,e