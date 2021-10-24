###############################################################################
# Figure S5
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
    # - ITNB_onset()
    # - S5_examples()

    
    
###############################################################################
# ITNB_onset
###############################################################################
    
def ITNB_onset(master_log, mice):
        
    temp_units1_totallicks=[]
    temp_mousename1=[]
    temp_KStest=[] 
    temp_unit_names=[]
    count=0
    temp_units1=[]
    temp_units2=[]           
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unit_names.append(unit_log['unit_name'].values[0])
                    temp_unit=[]
                    temp_unit_licks=[]
                    #print(unit)
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        temp_intrial=[] #to store spikes
                        temp_intrial_licks=[] # to store licks
                        for lick in trial_log['NR_Right_licks_aligned_spikes_NonBouts']:
                            temp_intrial.append(lick) 
                            temp_pre.append(np.sum((lick>-2) & (lick<-1)))
                            temp_post.append(np.sum((lick>-0.2) & (lick<0.8)))
                        temp_intrial_licks.append(len(trial_log['non_rewarded_Right_licksNonBouts'])) #get the numbe of right non-rewarded licks, weill be used to normalize
                        temp_unit.append(temp_intrial)
                        temp_unit_licks.append(temp_intrial_licks)
                    K,p=stats.ks_2samp(temp_pre,temp_post)
                    temp_KStest.append(p)
                    unit_spikes=[x for l in temp_unit for x in l]    #this listcpmprehension helps unpack the arrays.
                    temp_units1.append(unit_spikes) #list, one entry per unit, each entry has multiple 1d arrays of spikes
                    temp_units1_totallicks.append(temp_unit_licks)#list, one entry per unit, each entry has multiple scalars with number of licks in trial
                    temp_mousename1.append(mouse)
                    count+=1
    
    temp_units2_totallicks=[]
    temp_mousename2=[]
    temp_KStest=[] 
    count=0
    temp_units2=[]           
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unit=[]
                    temp_unit_licks=[]
                    #print(unit)
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        temp_intrial=[] #to store spikes
                        temp_intrial_licks=[] # to store licks
                        for lick in trial_log['NR_Left_licks_aligned_spikes_NonBouts']:
                            temp_intrial.append(lick) 
                            temp_pre.append(np.sum((lick>-2) & (lick<-1)))
                            temp_post.append(np.sum((lick>-0.2) & (lick<0.8)))
                        temp_intrial_licks.append(len(trial_log['non_rewarded_Left_licksNonBouts'])) #get the numbe of right non-rewarded licks, weill be used to normalize
                        temp_unit.append(temp_intrial)
                        temp_unit_licks.append(temp_intrial_licks)
                    K,p=stats.ks_2samp(temp_pre,temp_post)
                    temp_KStest.append(p)
                    unit_spikes=[x for l in temp_unit for x in l]    #this listcpmprehension helps unpack the arrays.
                    temp_units2.append(unit_spikes) #list, one entry per unit, each entry has multiple 1d arrays of spikes
                    temp_units2_totallicks.append(temp_unit_licks)#list, one entry per unit, each entry has multiple scalars with number of licks in trial
                    temp_mousename2.append(mouse)
                    count+=1
    
    
    peak_Right=[]
    peak_Right_norm=[]
    testR=[]
    for each in np.arange(len(temp_units1)):
        Spikes=temp_units1[each]
        Spikes_flat=[x for l in Spikes for x in l] #all 'non_rewarded_Right_licks' for unit
        Spikes_superflat=np.hstack(Spikes_flat)
        NumLicks=sum(np.asarray(temp_units1_totallicks[each])) #total number of licks that generated the spikes
        data,x=np.histogram(Spikes_superflat,bins=160, range=(-2,2) )
        data=data/NumLicks #Normalize
        testR.append(data)
        #determine whether min or max
        if np.mean(data[140:150])>np.mean(data[:20]):
            peak_Right.append(np.sum(data[np.argmax(data)-4:np.argmax(data)+4]) ) #This would fail if peak is close to the edges
            peak_Right_norm.append(np.sum(data[np.argmax(data)-4:np.argmax(data)+4]) / np.sum(data[:8])) 
        else:
            peak_Right.append(np.sum(data[np.argmin(data)-4:np.argmin(data)+4]) ) #This would fail if peak is close to the edges
            peak_Right_norm.append(np.sum(data[np.argmin(data)-4:np.argmin(data)+4]) / np.sum(data[:8])) 
        
    peak_Left=[]
    peak_Left_norm=[]
    testL=[]
    for each in np.arange(len(temp_units2)):
        Spikes=temp_units2[each]
        Spikes_flat=[x for l in Spikes for x in l] #all 'non_rewarded_Right_licks' for unit
        Spikes_superflat=np.hstack(Spikes_flat)
        NumLicks=sum(np.asarray(temp_units2_totallicks[each])) #total number of licks that generated the spikes
        data,x=np.histogram(Spikes_superflat,bins=160, range=(-2,2) )
        data=data/NumLicks #Normalize
        testL.append(data)
        #determine whether min or max
        if np.mean(data[140:150])>np.mean(data[:20]):
            peak_Left.append(np.sum(data[np.argmax(data)-4:np.argmax(data)+4]) ) #This would fail if peak is close to the edges
            peak_Left_norm.append(np.sum(data[np.argmax(data)-4:np.argmax(data)+4]) / np.sum(data[:8])) 
        else:
            peak_Left.append(np.sum(data[np.argmin(data)-4:np.argmin(data)+4]) ) #This would fail if peak is close to the edges
            peak_Left_norm.append(np.sum(data[np.argmin(data)-4:np.argmin(data)+4]) / np.sum(data[:8]))

    
    def f(l):
            for i in range(len(l)-4):
                if all(l[i:i+4]):
                    return True, i,
            return False, np.nan 
        
    x=[]
    y=[]
    for i,unit in enumerate(testR): 
        BaselineSomHit=np.mean(unit[:25])
        BaselineVisHit=np.mean(testL[i][:25])
        SDBaselineSomHit=np.std(unit[:25])
        SDBaselineVisHit=np.std(testL[i][:25])
        ThresholdSomHit= BaselineSomHit+3* SDBaselineSomHit
        ThresholdVisHit= BaselineVisHit+3* SDBaselineVisHit
        event, SomHitOnset = f(unit>ThresholdSomHit)
        event, VisHitOnset = f(testL[i]>ThresholdVisHit)
        x.append((SomHitOnset-80)*25) #time in milliseconds (neg is pre-lick)
        y.append((VisHitOnset-80)*25) #time in milliseconds
    
    Right_onset=x
    Left_onset=y    
    fig,ax = plt.subplots(1,2,figsize=(12,6))
    ax[0].set_aspect(aspect=1)
    #ax[0].fill([-1000,-1000,2000],[-1000,2000,2000], 'cornflowerblue', alpha=0.3)
    #ax[0].fill([-1000,2000,2000],[-1000,-1000,2000], 'gray', alpha=0.2)
    ax[0].plot(np.arange(-2000,2000), np.arange(-2000,2000), color='gray')   
    count=0
    for a,b in zip(x,y):
    #    if temp_unit_names[count] in TouchBlock_list:           
    #        ax[0].scatter(a,b, s=10, c='cornflowerblue')
    #    elif temp_unit_names[count] in VisualBlock_list:           
    #        ax[0].scatter(a,b, s=10, c='mediumseagreen')
    #    else:
            ax[0].scatter(a,b, s=10, c='k', alpha=0.5, linewidths=0)
        #count+=1
    Index= np.isnan(x) | np.isnan(y)
    Right_onset=[X for i,X in enumerate(x) if not Index[i]]
    Left_onset=  [X for i,X in enumerate(y) if not Index[i]]
    s,p=stats.wilcoxon( Right_onset, Left_onset)
    ax[0].set_title('ITNB-ActvityOnset', size=8)
    ax[0].set_ylabel('ActivityOnset \n LeftLicks ITNB (s)', size=20)
    ax[0].set_xlabel('ActivityOnset \n RightLicks ITNB (s)', size=20)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_xlim(-1200,2000)
    ax[0].set_ylim(-1200,2000)
    ax[0].set_xticks([-1000, 0, 1000, 2000])
    ax[0].set_xticklabels(['-1','0','1','2'], size=16)
    ax[0].set_yticks([-1000, 0, 1000, 2000])
    ax[0].set_yticklabels(['-1','0','1','2'], size=16)
    
    ax[1].set_aspect(aspect=3000)
    i = [i for i in x if ~np.isnan(i)]
    Cum=sp.stats.cumfreq(i, numbins=100)
    M=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    ax[1].plot(M,Cum.cumcount/np.size(i), color='k', label='RightLick')
    i = [i for i in y if ~np.isnan(i)]
    Cum=sp.stats.cumfreq(i, numbins=100)
    M=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    ax[1].plot(M,Cum.cumcount/np.size(i), color='grey', label='LeftLick')
    ax[1].set_xlim(-1200,2000)
    ax[1].set_ylim(0,1)
    ax[1].legend()
    ax[1].set_ylabel('CumFreq', size=20)          
    ax[1].set_title(' wilcoxon-p:'+str(p), size=8)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].set_xticks([-1000, 0, 1000, 2000])
    ax[1].set_xticklabels(['-1','0','1','2'], size=16)
    ax[1].set_xlabel('ActivityOnset-ITNB (s)', size=20)
    
    mean_right=np.mean(Right_onset)
    sem_right=np.std(Right_onset)/np.sqrt(len(Right_onset))
    mean_left=np.mean(Left_onset)
    sem_left=np.std(Left_onset)/np.sqrt(len(Left_onset))
    
    return mean_right, sem_right, mean_left, sem_left, temp_unit_names, temp_units1_totallicks, temp_units2_totallicks, temp_units1, temp_units2


###############################################################################
# S5_examples
###############################################################################
    
def S5_examples(master_log, mice, temp_units1, temp_units2, temp_unit_names,  temp_units1_totallicks, temp_units2_totallicks, directory): #if we do that, need to hand pick new examples
    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31','Claustrum32','Claustrum37']
    REVERSED=['Claustrum18','Claustrum23','Claustrum25']
    for example in [18,26,27,32,99,104,143,147,148,151,158,161,315,401,406,414,416,465]:
        fig,[ax1,ax2]=plt.subplots(2,1, figsize=(8,4))
        plt.subplots_adjust(wspace=0.2, hspace=0)
        plt.sca(ax1)
        Spikes1=temp_units1[example]
        for i in np.arange(len(Spikes1)):
            plt.scatter(Spikes1[i], np.zeros_like(Spikes1[i])+i, c='k', s=0.2)
        Spikes2=temp_units2[example]
        for i in np.arange(len(Spikes2)):
            plt.scatter(Spikes2[i], np.zeros_like(Spikes2[i])+i+len(Spikes1), c='grey', s=0.2)
        plt.xlim(-2,2)
        plt.vlines(0, 0,len(Spikes1)+len(Spikes2), color='r')
        plt.ylabel('Licks', size=16)
        plt.yticks([0,100,200],['0','100','200'], size=14)
        plt.xticks([],[])
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        
        Spikes_flat1=[x for l in Spikes1 for x in l] #all 'non_rewarded_Right_licks' for unit
        Spikes_flat2=[x for l in Spikes2 for x in l] #all 'non_rewarded_Right_licks' for unit
        
        Spikes_superflat1=np.hstack(Spikes_flat1)
        Spikes_superflat2=np.hstack(Spikes_flat2)
        
        NumLicks1=sum(np.asarray(temp_units1_totallicks[example])) #total number of licks that generated the spikes
        NumLicks2=sum(np.asarray(temp_units2_totallicks[example])) #total number of licks that generated the spikes
        
        data1,x=np.histogram(Spikes_superflat1,bins=100, range=(-2,2) )
        data1=data1/NumLicks1 #Normalize
        data1=data1/np.mean(data1[:25])
        data2,x=np.histogram(Spikes_superflat2,bins=100, range=(-2,2) )
        data2=data2/NumLicks2 #Normalize
        data2=data2/np.mean(data2[:25])
        
        plt.sca(ax2)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        plt.plot(np.arange(100), data1, color='k')
        plt.plot(np.arange(100), data2, color='grey')
        plt.xticks([0,25,50,75,100],['-2','-1','0','1','2'], size=14)
        plt.vlines(50, 0,6, color='r')
        plt.ylabel('Activity\n(Norm)', size=16)
        plt.xlabel('Time (seconds)', size=16)
        plt.xlim(0,100)
        plt.yticks([0,2,4,6],['0','2','4', '6'], size=14)
        plt.savefig(directory+temp_unit_names[example]+'_ITNB.png', bbox_inches='tight')
        plt.close()
        #########################################
        fig,[ax3,ax4]=plt.subplots(2,1, figsize=(8,4))
        plt.subplots_adjust(wspace=0.2, hspace=0)
        
        unit_log=master_log.loc[np.equal(master_log['unit_name'], temp_unit_names[example])]
        if unit_log['mouse_name'].values[0] in ORIGINAL:
            color1='#A1007D'
            color2='#00D88A'
        elif unit_log['mouse_name'].values[0] in REVERSED:
            color1='#00D88A'
            color2='#A1007D'
        SomHit_spikes=[]
        SomHit_licks=[]
        for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
            trial_log=unit_log.loc[trial,:]
            SomHit_spikes.append(trial_log['LickALigned_spike_times'])
            SomHit_licks.append(trial_log['FirstLick'][0][0])
            
        VisHit_spikes=[]
        VisHit_licks=[]
        for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
            trial_log=unit_log.loc[trial,:]
            VisHit_spikes.append(trial_log['LickALigned_spike_times'])
            VisHit_licks.append(trial_log['FirstLick'][0][0])
            
        
        plt.sca(ax3)
        counter=0
        SomHitindex=np.argsort(SomHit_licks)
        for i in SomHitindex:
            spikes=SomHit_spikes[i]
            ax3.scatter(spikes, np.zeros_like(spikes)+counter, c=color1, s=0.5) 
            ax3.plot(0-SomHit_licks[i],counter, color='grey',marker='.', markersize=2)
            counter+=1
        
        ax3.set_xlim(-2,2)
        ax3.set_xticks([-1,0,1])
        ax3.set_xticklabels([])
        ax3.vlines(0,1,counter-1, color='r')
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.set_ylabel('Trials', size=16)
        
        VisHitindex=np.argsort(VisHit_licks)   
        for i in VisHitindex:
            spikes=VisHit_spikes[i]
            ax3.scatter(spikes, np.zeros_like(spikes)+counter, c=color2, s=0.5) 
            ax3.plot(0-VisHit_licks[i],counter, color='grey',marker='.', markersize=2)
            counter+=1
        
        ax3.set_xlim(-2,2)
        ax3.set_xticks([-1,0,1])
        ax3.set_xticklabels([-1,0,1])
        ax3.vlines(0,1,counter-1, color='r')
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        
        
        data1,x=np.histogram([x for l in SomHit_spikes for x in l],bins=100, range=(-2,2) )
        data1=data1/NumLicks1 #Normalize
        data1=data1/np.mean(data1[:25])
        data2,x=np.histogram([x for l in VisHit_spikes for x in l],bins=100, range=(-2,2) )
        data2=data2/NumLicks2 #Normalize
        data2=data2/np.mean(data2[:25])
        
        plt.sca(ax4)
        ax4.spines['right'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        plt.plot(np.arange(100), data1, color=color1)
        plt.plot(np.arange(100), data2, color=color2)
        plt.xticks([0,25,50,75,100],['-2','-1','0', '1','2'], size=14)
        plt.vlines(50, 0,6, color='r')
        plt.ylabel('Activity\n(Norm)', size=16)
        plt.xlabel('Time (s)', size=16)
        plt.xlim(0,100)
        plt.yticks([0,2,4,6],['0','2','4', '6'], size=14)
        plt.savefig(directory+temp_unit_names[example]+'_Task.png', bbox_inches='tight')
        plt.close()

    return