###############################################################################
# Figure S2
###############################################################################

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
from Helper_functions import MidpointNormalize

# FUNCTIONS: 
    # - example_claustrum_lickaligned
    # - Heatmap_AllTrialTypes_lickaligned
    # - S1_CRresponse
    # - Stim_Lick_aligned
    
##############################################################################
# example_claustrum
##############################################################################

def example_claustrum_lickaligned(unit_log, yaxis_range):
       
    SomHit_spikes=[]
    SomHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
        trial_log=unit_log.loc[trial,:]
        SomHit_spikes.append(trial_log['LickALigned_spike_times'])
        SomHit_licks.append(trial_log['FirstLick'][0][0])
    
    SomFA_spikes=[]
    SomFA_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
        trial_log=unit_log.loc[trial,:]
        SomFA_spikes.append(trial_log['LickALigned_spike_times'])
        SomFA_licks.append(trial_log['FirstLick'][0][0])
    
    VisHit_spikes=[]
    VisHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
        trial_log=unit_log.loc[trial,:]
        VisHit_spikes.append(trial_log['LickALigned_spike_times'])
        VisHit_licks.append(trial_log['FirstLick'][0][0])
    
    VisFA_spikes=[]
    VisFA_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisFA'].index:
        trial_log=unit_log.loc[trial,:]
        VisFA_spikes.append(trial_log['LickALigned_spike_times'])
        VisFA_licks.append(trial_log['FirstLick'][0][0])
    
    #--------------------#
        #FIG
    #--------------------#
    
    fig, ax=plt.subplots(4,1, figsize=(6,16))
    plt.subplots_adjust(wspace=0.2, hspace=0.05)
    ax=ax.flatten()
    
    counter=0
    SomHitindex=np.argsort(SomHit_licks)
    SomFAindex=np.argsort(SomFA_licks)
    
    for i in SomFAindex:
        spikes=SomFA_spikes[i]
        ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='tomato', s=0.5) 
        ax[0].plot(-SomFA_licks[i],counter, 'k.', markersize=2)
        counter+=1
        
    
    for i in SomHitindex:
        if SomHit_licks[i]>0.1:
            spikes=SomHit_spikes[i]
            ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
            ax[0].plot(-SomHit_licks[i],counter, 'k.', markersize=2)
            counter+=1
    
    ax[0].set_xlim(-2,2)
    ax[0].set_xticks([-2,-1, 0,1,2])
    ax[0].set_xticklabels([])
    ax[0].vlines(0,1,counter-1)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_ylabel('Trials', size=18)
    
    counter=0
    VisHitindex=np.argsort(VisHit_licks)
    VisFAindex=np.argsort(VisFA_licks)
    
    for i in VisHitindex:
        if VisHit_licks[i]>0.1:
            spikes=VisHit_spikes[i]
            ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='tomato', s=0.5) 
            ax[2].plot(-VisHit_licks[i],counter, 'k.', markersize=2)
            counter+=1
    
    for i in VisFAindex:
        spikes=VisFA_spikes[i]
        ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
        ax[2].plot(-VisFA_licks[i],counter, 'k.', markersize=2)
        counter+=1
    ax[2].set_xlim(-2,2)
    ax[2].set_xticks([-2,-1,0,1,2])
    ax[2].set_xticklabels([])
    ax[2].vlines(0,1,counter-1)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].set_ylabel('Trials', size=18)
    
    #mean
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
    temp_units2=np.mean(temp_trial, axis=0)
    
    
    
    #temp_units1_sma=pd.DataFrame(np.transpose(temp_units1)).rolling(window=10, center=True).mean()
    #temp_units2_sma=pd.DataFrame(np.transpose(temp_units2)).rolling(window=10, center=True).mean()
    #temp_units3_sma=pd.DataFrame(np.transpose(temp_units3)).rolling(window=10, center=True).mean()
    #temp_units4_sma=pd.DataFrame(np.transpose(temp_units4)).rolling(window=10, center=True).mean()
    
    ax[1].plot(np.arange(79), (temp_units2[0:-1:2] + temp_units2[1::2]) *20, color='tomato', alpha=0.7)
    ax[1].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='slateblue', alpha=0.7)
    
    ax[1].set_xticks([19.5, 39.5, 59.5])
    ax[1].set_xticklabels([])
       # ax[1].set_xlabel('Time (s)', size=18)
    ax[1].set_ylabel('Activity (Hz)', size=18)
    ax[1].vlines(39.5,0,yaxis_range[1])
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].set_ylim(0,yaxis_range[1])
    ax[1].set_xlim(-0.5,79.5)
    #
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisFA'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
    temp_units2=np.mean(temp_trial, axis=0)
    
    
    ax[3].plot(np.arange(79), (temp_units2[0:-1:2] + temp_units2[1::2]) *20, color='slateblue', alpha=0.7)
    ax[3].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='tomato', alpha=0.7)
    
    ax[3].set_xticks([-0.5,19.5, 39.5, 59.5, 79.5])
    ax[3].set_xticklabels(['-2','-1','0','1','2'])
    ax[3].vlines(39.5,0,yaxis_range[1])
    ax[3].set_xlabel('Time (s)', size=18)
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    ax[3].set_ylim(0,yaxis_range[1])
    ax[3].set_xlim(-0.5,79.5)
    ax[3].set_ylabel('Activity (Hz)', size=18)
    
    fig.suptitle(unit_log['unit_name'].values[0])
    return


###############################################################################
# Heatmap_AllTrialTypes_lickaligned
###############################################################################
    
def Heatmap_AllTrialTypes_lickaligned(mice, master_log):
    from Helper_functions import MidpointNormalize
    plt.style.use('default')
    fig, ax=plt.subplots(2,2,figsize=(6,8))   
    ax=ax.flatten()
    trialtypes=['SomHit','SomFA', 'VisFA', 'VisHit']
    plot_pos=[1,0,2,3]
    for i,trialtype in zip(plot_pos,trialtypes):
        temp_units1=[]
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    norm=np.mean(np.mean(unit_log['-2to2sec_25msecbins_LickAligned'],0)[:39])
                    temp_trial=[]
                    for trial in unit_log[unit_log['Stim/Block/Response']==trialtype].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
                    if norm:
                        temp_units1.append(np.mean(temp_trial, axis=0)/norm)
                    else:
                        temp_units1.append(np.mean(temp_trial, axis=0))
        if trialtype=='SomHit':
            order=np.argsort([np.sum(x[40:80]) for x in temp_units1])
        sorted_temp_units1=[temp_units1[row] for row in order[::-1]]
        
        sorted_temp_units1=[x if np.nansum(x) else np.zeros(159) for x in sorted_temp_units1]
        
        SomHit=ax[i-1].imshow(sorted_temp_units1,norm=MidpointNormalize(midpoint=1.,vmin=0, vmax=4), cmap='bwr')
        
        ax[i-1].vlines(79.5,0,len(temp_units1)-1, color='k', linewidth=2)
        ax[i-1].set_title('SomHit')
        ax[i-1].set_xticks([39.5,79.5,119.5])
        ax[i-1].set_xticklabels([])
        ax[i-1].set_yticks([])
        ax[i-1].set_yticklabels([])
        ax[i-1].set_title(trialtype)
    
        if (trialtype=='SomHit') | (trialtype=='VisHit'):
          ax[i-1].set_yticks([0,50,100,150,200])
          ax[i-1].set_yticklabels(['0','50','100','150','200'])
          ax[i-1].set_ylabel('Normalized activity')
          
        if (trialtype=='VisHit') | (trialtype=='SomFA') :
          ax[i-1].set_xticks([-0.05,39.5, 79.5, 119.5, 159.5])
          ax[i-1].set_xticklabels(['-2','-1','0', '1','2'])
          ax[i-1].set_xlabel('Time (seconds)')
          
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.1)
    # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
    cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(SomHit, cax=cb_ax)
    return


###############################################################################
# S1_CRresponse
###############################################################################

def S1_CRresponse(S1master_log, mice):
    
    # Are CR responses in S1 suppressed compared to Hit
    
    mice=['EF0074', 'EF0076','EF0077','EF0079']#,'EF0081','EF0085', 'EF0083', 'EF0084','EF0088','EF0089',]
    
    plt.style.use('default')
    fig, ax=plt.subplots(1,1,figsize=(12,8))   
    trialtypes=['SomHit','VisCR']
    colors=['slateblue','cornflowerblue']
    Hit_data=[]
    CR_data=[]
    for i in np.arange(len(trialtypes)):
        temp_units1=[]
        for mouse in mice:
            mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
                    temp_trial=[]
                    for trial in unit_log[unit_log['Stim/Block/Response']==trialtypes[i]].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned'])
                    
                    if trialtypes[i]=='SomHit':
                        Hit_data.append(np.mean(temp_trial, axis=0)[41])
                    elif trialtypes[i]=='VisCR':
                        CR_data.append(np.mean(temp_trial, axis=0)[41])
                    if norm:
                        temp_units1.append(np.mean(temp_trial, axis=0)/norm)
                    else:
                        temp_units1.append(np.mean(temp_trial, axis=0))
        SEM1=sp.stats.sem(temp_units1, axis=0)
        SomHit, =ax.plot(np.arange(159), np.mean(temp_units1, axis=0), color=colors[i])
        ax.fill_between(np.arange(159),np.mean(temp_units1, axis=0)+SEM1, np.mean(temp_units1, axis=0)-SEM1, color=colors[i], alpha=0.1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.vlines(39,0.8,2)         
    plt.xticks([39,41,43,45,47,49], ['0','50','100','150','200','250'], size=14)
    plt.yticks([1,1.5,2], size=14)
    plt.xlabel('Time (ms))', size=16)
    plt.ylabel('Normalized activity', size=16)
    plt.xlim(38,50)   
    plt.title('SomHit vs VisCR S1, n='+str(len(temp_units1)))
    
    fig, ax=plt.subplots(1,1,figsize=(8,8))
    s,p=sp.stats.wilcoxon(Hit_data, CR_data)
    ax.scatter( [np.log(x*40) for x in CR_data], [np.log(x*40) for x in Hit_data], c='cornflowerblue', alpha=0.8, linewidths=0, s=30)
    ax.plot([-1,np.log(150)], [-1,np.log(150)], color='k', linestyle='dotted')
    meanHIT=np.log(np.mean(Hit_data)*40)
    semHIT=np.log(np.std(Hit_data)*40/np.sqrt(len(Hit_data)))
    meanCR=np.log(np.mean(CR_data)*40)
    semCR=np.log(np.std(CR_data)*40/np.sqrt(len(CR_data)))
    plt.hlines(meanHIT, meanCR-semCR, meanCR+semCR)
    plt.vlines(meanCR, meanHIT-semHIT, meanHIT+semHIT)
    plt.xlim(-1,np.log(150))
    plt.ylim(-1,np.log(150))

    plt.xticks([np.log(x) if x!=0 else 0 for x in [0,5,10,20,40,80]],['0','5','10','20','40', '80'],  size=16)
    plt.yticks([np.log(x) if x!=0 else 0 for x in [0,5,10,20,40, 80]],['0','5','10','20','40', '80'],  size=16)
    plt.ylabel('Attended Touch Response (Hz)', size=20)
    plt.xlabel('Ignored Touch Response (Hz)', size=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.title('SomHit vs VisCR S1 - wilcoxon p:'+str(p))
     
    #Numbers
    Mean_Hit=np.mean([x*40 for x in Hit_data])
    SEM_Hit=np.std([x*40 for x in Hit_data])/np.sqrt(len(Hit_data))
    Mean_CR=np.mean([x*40 for x in CR_data])
    SEM_CR=np.std([x*40 for x in CR_data])/np.sqrt(len(CR_data))

    return Mean_Hit, SEM_Hit, Mean_CR, SEM_CR

###############################################################################
# Stim_Lick_aligned
###############################################################################

def Stim_Lick_aligned(master_log, trial_type, example):
    unit_log=master_log[master_log['unit_name']==example]
    yaxis_range=[0,50]
     
    SomHit_spikes=[]
    SomHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        SomHit_spikes.append(trial_log['StimALigned_spike_times'])
        SomHit_licks.append(trial_log['FirstLick'][0][0])
        
       
    
    fig, ax=plt.subplots(2,2, figsize=(6,6))
    plt.subplots_adjust(wspace=0.2, hspace=0.05)
    ax=ax.flatten()
    
    counter=0
    SomHitindex=np.argsort(SomHit_licks)
       
    for i in SomHitindex:
        if SomHit_licks[i]>0.1:
            spikes=SomHit_spikes[i]
            ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
            ax[0].plot(SomHit_licks[i],counter, 'k.', markersize=2)
            counter+=1
    
    ax[0].set_xlim(-1,3)
    ax[0].set_xticks([-1,0,1,2,3])
    ax[0].set_xticklabels([])
    ax[0].vlines(0,1,counter-1)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_ylabel('Trials', size=18)
    
       
    #mean
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
      
    ax[2].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='slateblue', alpha=0.7)
    
    ax[2].set_xticks([0, 19.5, 39.5, 59.5, 79.5])
    ax[2].set_xticklabels(['-1','0','1','2', '3'])
    ax[2].set_xlabel('Time (s)', size=18)
    ax[2].set_ylabel('Activity (Hz)', size=18)
    ax[2].vlines(19.5,0,yaxis_range[1])
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
   # ax[2].set_ylim(0,yaxis_range[1])
    ax[2].set_xlim(-0.5,79.5)
       
    
    SomHit_spikes=[]
    SomHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        SomHit_spikes.append(trial_log['LickALigned_spike_times'])
        SomHit_licks.append(trial_log['FirstLick'][0][0])
        
       
    
    
    counter=0
    SomHitindex=np.argsort(SomHit_licks)
       
    for i in SomHitindex:
        if SomHit_licks[i]>0.1:
            spikes=SomHit_spikes[i]
            ax[1].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
            ax[1].plot(SomHit_licks[i]*(-1),counter, 'k.', markersize=2)
            counter+=1
    
    ax[1].set_xlim(-2,2)
    ax[1].set_xticks([-2,-1,0,1,2])
    ax[1].set_xticklabels([])
    ax[1].vlines(0,1,counter-1)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    
       
    #mean
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
      
    ax[3].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='slateblue', alpha=0.7)
    
    ax[3].set_xticks([0, 19.5, 39.5, 59.5, 79.5])
    ax[3].set_xticklabels(['-2','-1','0','1','2'])
    ax[3].set_xlabel('Time (s)', size=18)
    ax[3].vlines(39.5,0,yaxis_range[1])
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
   # ax[3].set_ylim(0,yaxis_range[1])
    ax[3].set_xlim(-0.5,79.5)
    
    fig.suptitle(unit_log['unit_name'].values[0])
#    plt.savefig('C:\\Users\\Brown Lab\\Desktop\\trash\\' + example +'.jpg')
#    plt.close('all')
    return