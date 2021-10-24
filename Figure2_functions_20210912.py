###############################################################################
# Figure 2
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
    # - example_claustrum
    # - Heatmap_AllTrialTypes
    # - Responsiveness
    
##############################################################################
# example_claustrum
##############################################################################

def example_claustrum(unit_log, yaxis_range):
    ###########
    #Example
    ###########
    
    
    norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
    
    SomHit_spikes=[]
    SomHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
        trial_log=unit_log.loc[trial,:]
        SomHit_spikes.append(trial_log['StimALigned_spike_times'])
        SomHit_licks.append(trial_log['FirstLick'][0][0])
        
    SomMiss_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomMiss'].index:
        trial_log=unit_log.loc[trial,:]
        SomMiss_spikes.append(trial_log['StimALigned_spike_times'])
    
    SomFA_spikes=[]
    SomFA_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
        trial_log=unit_log.loc[trial,:]
        SomFA_spikes.append(trial_log['StimALigned_spike_times'])
        SomFA_licks.append(trial_log['FirstLick'][0][0])
    
    SomCR_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomCR'].index:
        trial_log=unit_log.loc[trial,:]
        SomCR_spikes.append(trial_log['StimALigned_spike_times'])
        
    VisHit_spikes=[]
    VisHit_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
        trial_log=unit_log.loc[trial,:]
        VisHit_spikes.append(trial_log['StimALigned_spike_times'])
        VisHit_licks.append(trial_log['FirstLick'][0][0])
        
    VisMiss_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisMiss'].index:
        trial_log=unit_log.loc[trial,:]
        VisMiss_spikes.append(trial_log['StimALigned_spike_times'])
    
    VisFA_spikes=[]
    VisFA_licks=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisFA'].index:
        trial_log=unit_log.loc[trial,:]
        VisFA_spikes.append(trial_log['StimALigned_spike_times'])
        VisFA_licks.append(trial_log['FirstLick'][0][0])
    
    VisCR_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisCR'].index:
        trial_log=unit_log.loc[trial,:]
        VisCR_spikes.append(trial_log['StimALigned_spike_times'])
        
    #--------------------#
        #FIG
    #--------------------#
    
    fig, ax=plt.subplots(4,1, figsize=(6,16))
    plt.subplots_adjust(wspace=0.2, hspace=0.05)
    ax=ax.flatten()
    
    counter=0
    SomHitindex=np.argsort(SomHit_licks)
    SomFAindex=np.argsort(SomFA_licks)
    for i in np.arange(len(SomCR_spikes)):
        spikes=SomCR_spikes[i]
        ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='orange', s=0.5) 
        counter+=1
    
    for i in SomFAindex:
        spikes=SomFA_spikes[i]
        ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='tomato', s=0.5) 
        ax[0].plot(SomFA_licks[i],counter, 'k.', markersize=2)
        counter+=1
        
    for i in np.arange(len(SomMiss_spikes)):
        spikes=SomMiss_spikes[i]
        ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='cornflowerblue', s=0.5) 
        counter+=1
    
    for i in SomHitindex:
        if SomHit_licks[i]>0.1:
            spikes=SomHit_spikes[i]
            ax[0].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
            ax[0].plot(SomHit_licks[i],counter, 'k.', markersize=2)
            counter+=1
    
    ax[0].set_xlim(-1,3)
    ax[0].set_xticks([0,1,2])
    ax[0].set_xticklabels([])
    ax[0].vlines(0,1,counter-1)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_ylabel('Trials', size=18)
    
    counter=0
    VisHitindex=np.argsort(VisHit_licks)
    VisFAindex=np.argsort(VisFA_licks)
    for i in np.arange(len(VisMiss_spikes)):
        spikes=VisMiss_spikes[i]
        ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='orange', s=0.5) 
        counter+=1
        
    for i in VisHitindex:
        if VisHit_licks[i]>0.1:
            spikes=VisHit_spikes[i]
            ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='tomato', s=0.5) 
            ax[2].plot(VisHit_licks[i],counter, 'k.', markersize=2)
            counter+=1
        
    for i in np.arange(len(VisCR_spikes)):
        spikes=VisCR_spikes[i]
        ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='cornflowerblue', s=0.5) 
        counter+=1
        
    for i in VisFAindex:
        spikes=VisFA_spikes[i]
        ax[2].scatter(spikes, np.zeros_like(spikes)+counter, c='slateblue', s=0.5) 
        ax[2].plot(VisFA_licks[i],counter, 'k.', markersize=2)
        counter+=1
    ax[2].set_xlim(-1,3)
    ax[2].set_xticks([0,1,2])
    ax[2].set_xticklabels([])
    ax[2].vlines(0,1,counter-1)
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[2].set_ylabel('Trials', size=18)
    
    #mean
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units2=np.mean(temp_trial, axis=0)
    
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomMiss'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units3=np.mean(temp_trial, axis=0)
    
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='SomCR'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units4=np.mean(temp_trial, axis=0)
    
    #temp_units1_sma=pd.DataFrame(np.transpose(temp_units1)).rolling(window=10, center=True).mean()
    #temp_units2_sma=pd.DataFrame(np.transpose(temp_units2)).rolling(window=10, center=True).mean()
    #temp_units3_sma=pd.DataFrame(np.transpose(temp_units3)).rolling(window=10, center=True).mean()
    #temp_units4_sma=pd.DataFrame(np.transpose(temp_units4)).rolling(window=10, center=True).mean()
    
    ax[1].plot(np.arange(79), (temp_units2[0:-1:2] + temp_units2[1::2]) *20, color='tomato', alpha=0.7)
    ax[1].plot(np.arange(79),(temp_units3[0:-1:2] + temp_units3[1::2]) *20, color='cornflowerblue', alpha=0.7)
    ax[1].plot(np.arange(79),(temp_units4[0:-1:2] + temp_units4[1::2]) *20, color='orange', alpha=0.7)
    ax[1].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='slateblue', alpha=0.7)
    
    ax[1].set_xticks([19.5, 39.5, 59.5])
    ax[1].set_xticklabels([])
   # ax[1].set_xlabel('Time (s)', size=18)
    ax[1].set_ylabel('Activity (Hz)', size=18)
    ax[1].vlines(19.5,0,yaxis_range[1])
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].set_ylim(0,yaxis_range[1])
    ax[1].set_xlim(-0.5,79.5)
    #
    temp_trial=[]
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units1=np.mean(temp_trial, axis=0)
        
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisFA'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units2=np.mean(temp_trial, axis=0)
    
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisMiss'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units3=np.mean(temp_trial, axis=0)
    
    temp_trial=[]    
    for trial in unit_log[unit_log['Stim/Block/Response']=='VisCR'].index:
        trial_log=unit_log.loc[trial,:]
        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    temp_units4=np.mean(temp_trial, axis=0)
    
    
    ax[3].plot(np.arange(79), (temp_units2[0:-1:2] + temp_units2[1::2]) *20, color='slateblue', alpha=0.7)
    ax[3].plot(np.arange(79),(temp_units3[0:-1:2] + temp_units3[1::2]) *20, color='orange', alpha=0.7)
    ax[3].plot(np.arange(79),(temp_units4[0:-1:2] + temp_units4[1::2]) *20, color='cornflowerblue', alpha=0.7)
    ax[3].plot(np.arange(79),(temp_units1[0:-1:2] + temp_units1[1::2]) *20, color='tomato', alpha=0.7)
    
    ax[3].set_xticks([-0.5,19.5, 39.5, 59.5, 79.5])
    ax[3].set_xticklabels(['-1','0','1','2','3'])
    ax[3].vlines(19.5,0,yaxis_range[1])
    ax[3].set_xlabel('Time (s)', size=18)
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    ax[3].set_ylim(0,yaxis_range[1])
    ax[3].set_xlim(-0.5,79.5)
    ax[3].set_ylabel('Activity (Hz)', size=18)
    
    fig.suptitle(unit_log['unit_name'].values[0])
    return

###############################################################################
# Heatmap_AllTrialTypes
###############################################################################
    
def Heatmap_AllTrialTypes(mice, master_log):
    #mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37'] 
    

    plt.style.use('default')
    fig, ax=plt.subplots(2,4,figsize=(12,8))   
    ax=ax.flatten()
    trialtypes=['SomHit','VisMiss','SomFA','SomMiss','SomCR', 'VisFA', 'VisHit','VisCR']
    plot_pos=[1,0,2,3,4,5,6,7]
    for i,trialtype in zip(plot_pos,trialtypes):
        temp_units1=[]
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
                    temp_trial=[]
                    for trial in unit_log[unit_log['Stim/Block/Response']==trialtype].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
                    if norm:
                        temp_units1.append(np.mean(temp_trial, axis=0)/norm)
                    else:
                        temp_units1.append(np.mean(temp_trial, axis=0))
        if trialtype=='SomHit':
            order=np.argsort([np.sum(x[40:80]) for x in temp_units1])
        sorted_temp_units1=[temp_units1[row] for row in order[::-1]]
        
        sorted_temp_units1=[x if np.nansum(x) else np.zeros(159) for x in sorted_temp_units1]
        
        SomHit=ax[i-1].imshow(sorted_temp_units1,norm=MidpointNormalize(midpoint=1.,vmin=0, vmax=4), cmap='bwr')
        
        ax[i-1].vlines(39.5,0,len(temp_units1)-1, color='k', linewidth=2)
        ax[i-1].set_title('SomHit')
        ax[i-1].set_xticks([39.5,79.5,119.5])
        ax[i-1].set_xticklabels([])
        ax[i-1].set_yticks([])
        ax[i-1].set_yticklabels([])
        ax[i-1].set_title(trialtype)
        
        if False in master_log['Category'].values:#added this condtion to only add the 200 tick marks for the full dataset, not the optoonly
            if (trialtype=='SomHit') | (trialtype=='VisFA'):
              ax[i-1].set_yticks([0,50,100,150,200])
              ax[i-1].set_yticklabels(['0','50','100','150','200'])
              ax[i-1].set_ylabel('Normalized activity')
          
        if (trialtype=='VisHit') | (trialtype=='VisFA') |  (trialtype=='VisMiss') | (trialtype=='VisCR'):
          ax[i-1].set_xticks([-0.05,39.5, 79.5, 119.5, 159.5])
          ax[i-1].set_xticklabels(['-1','0','1', '2','3'])
          ax[i-1].set_xlabel('Time (seconds)')
          
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.1)
    # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
    cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(SomHit, cax=cb_ax)
    return

###############################################################################
# Responsiveness
###############################################################################

def Responsiveness(trial_type, master_log, mice):
    
    if trial_type=='SomCR':
        color='orange'
    elif trial_type== 'VisCR':
        color='cornflowerblue'
    elif trial_type=='SomHit':
        color='slateblue'
    elif trial_type=='VisHit':
        color='tomato'
        
    unit_names=[]
    p_values=[]
    Direction=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        print(mouse)
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            print(day)
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                #print(unit_log['unit_name'].values[0])
                temp_PRE=[]
                temp_POST=[]
                for trial in unit_log[(unit_log['Stim/Block/Response']==trial_type) ].index:
                    trial_log=unit_log.loc[trial,:]
                    temp_PRE.append(np.sum((trial_log['StimALigned_spike_times']<0)&(trial_log['StimALigned_spike_times']>-1)))
                    temp_POST.append(np.sum((trial_log['StimALigned_spike_times']>0)&(trial_log['StimALigned_spike_times']<1)))
                s,p=sp.stats.ks_2samp(temp_PRE, temp_POST)
                Direction.append(np.mean(temp_POST)-np.mean(temp_PRE))
                unit_names.append(trial_log['unit_name'])
                p_values.append(p)
                

    Touch_responsive_units=[]
    for p,name,d in zip(p_values, unit_names,Direction):
        if p<0.05:
            Touch_responsive_units.append(name)
    
      
#    #loop used to look at all examples        
#     for j in np.arange( 400, 450):
#         unit_log=master_log.loc[np.equal(master_log['unit_name'], unit_names[j])]
#         norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
#        
##         Som_spikes=[]
##         for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
##             trial_log=unit_log.loc[trial,:]
##             Som_spikes.append(trial_log['StimALigned_spike_times'])
#      
#            
##         fig, ax=plt.subplots(1,1, figsize=(8,8))
##         for i,spikes in enumerate(Som_spikes):
##             ax.scatter(spikes, np.zeros_like(spikes)+i, c=color, s=0.5) 
##         ax.vlines(0,0,len(Som_spikes))
##         ax.set_xlim(-1,3)
##         ax.set_title(Touch_responsive_units[j])
#        
#         Som_spikes=[]
#         for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
#             trial_log=unit_log.loc[trial,:]
#             Som_spikes.append(trial_log['-1to3sec_25msecbins_StimAligned'])
#            
#         fig, ax=plt.subplots(1,1, figsize=(8,8))
#         ax.bar(np.arange(159), np.mean(Som_spikes, axis=0))
#         ax.vlines(40,0,2*np.max(np.mean(Som_spikes, axis=0)))
#         ax.set_title(unit_names[j])
#    
    #select and plot example
    if trial_type=='SomCR':
        unit_name='Claustrum37_20191123_TT8clst2'
    elif trial_type== 'VisCR':
        unit_name='Cl5_05-31-17_TT4clst2'
    elif trial_type=='SomHit':
        unit_name='Cl4_05-25-17_TT4clst3'
    elif trial_type=='VisHit':
        unit_name='Cl4_05-25-17_TT4clst3'
    
    unit_log=master_log.loc[np.equal(master_log['unit_name'], unit_name)]
    
    Som_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        Som_spikes.append(trial_log['StimALigned_spike_times'])
    
    fig, [ax1, ax2]=plt.subplots(2,1, figsize=(8,8))
    plt.subplots_adjust(wspace=0.2, hspace=0)
    for i,spikes in enumerate(Som_spikes):
        ax1.scatter(spikes, np.zeros_like(spikes)+i, c=color, s=3) 
    ax1.vlines(0,0,len(Som_spikes))
    ax1.set_xlim(-1,3)
    ax1.set_xticklabels([])
    ax1.set_title(unit_name)
    
    ax1.set_ylabel('trials', size=14)
    
    Som_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        Som_spikes.append(trial_log['-1to3sec_25msecbins_StimAligned'])
        
    ax2.bar(np.arange(159), np.mean(Som_spikes, axis=0)*40, color=color)
    ax2.vlines(39.5,0,2*40*np.max(np.mean(Som_spikes, axis=0)))
    ax2.set_xlim(-0.5, 159)
    ax2.set_ylim(0,45)
    ax2.set_xticks([-0.5, 39.5, 79.5, 119.5, 159.5])
    ax2.set_xticklabels(['-1','0','1','2','3'])
    ax2.set_xlabel('Time(seconds)', size=14)
    ax2.set_ylabel('spikes/sec', size=14)
    
    #select another neuron for negative example
    if trial_type=='SomCR':
        unit_name='Cl5_05-31-17_TT4clst6'
        unit_name='Claustrum23_20190731_TT6clst3'
    elif trial_type== 'VisCR':
        unit_name='Claustrum37_20191124_TT3clst1'
    elif trial_type=='SomHit':
        unit_name='Claustrum18_20190418_TT7clst1'
    elif trial_type=='VisHit':
        unit_name='Cl4_05-19-17_TT7clst1'
    
    unit_log=master_log.loc[np.equal(master_log['unit_name'], unit_name)]
    
    Som_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        Som_spikes.append(trial_log['StimALigned_spike_times'])
    
    fig, [ax1, ax2]=plt.subplots(2,1, figsize=(8,8))
    plt.subplots_adjust(wspace=0.2, hspace=0)
    for i,spikes in enumerate(Som_spikes):
        ax1.scatter(spikes, np.zeros_like(spikes)+i, c='grey', s=3) 
    ax1.vlines(0,0,len(Som_spikes))
    ax1.set_xlim(-1,3)
    ax1.set_xticklabels([])
    ax1.set_title(unit_name)
    
    ax1.set_ylabel('trials', size=14)
    
    Som_spikes=[]
    for trial in unit_log[unit_log['Stim/Block/Response']==trial_type].index:
        trial_log=unit_log.loc[trial,:]
        Som_spikes.append(trial_log['-1to3sec_25msecbins_StimAligned'])
        
    ax2.bar(np.arange(159), np.mean(Som_spikes, axis=0)*40, color='grey')
    ax2.vlines(39.5,0,45)
    ax2.set_xlim(-0.5, 159)
    ax2.set_ylim(0,45)
    ax2.set_xticks([-0.5,39.5, 79.5, 119.5, 159.5])
    ax2.set_xticklabels(['-1','0','1','2','3'])
    ax2.set_xlabel('Time(seconds)', size=14)
    ax2.set_ylabel('spikes/sec', size=14)
    
    ### Piechart of responsive units
    plt.figure()
    plt.pie([len(Touch_responsive_units), len(np.unique(master_log['unit_name']))-len(Touch_responsive_units)], 
            colors=[color, 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[len(Touch_responsive_units), len(np.unique(master_log['unit_name']))-len(Touch_responsive_units)],
            labeldistance=0.65 )
    plt.title('RESPONSIVE UNITS')
    

    
    Number_Resp= len(Touch_responsive_units)
    Number_Total=len(np.unique(master_log['unit_name']))
    
    return Number_Resp, Number_Total