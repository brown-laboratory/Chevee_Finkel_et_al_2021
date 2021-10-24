###############################################################################
# Figure 4
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
    # - AUC_All_LickDirection()
    # - ITNB_scatter()
    # - example_ITBNandAUC()
    # - R2_correlation()
    # - Loglikelihood_decoder()
    
    
    



###############################################################################
# Split by Lick Direction
###############################################################################

def AUC_All_LickDirection(mice, master_log, BinNumber, example1, example2):
    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
    REVERSED=['Claustrum18','Claustrum23','Claustrum25'] 
    
    temp_unitsA=[]
    for mouse in mice:
        if mouse in ORIGINAL:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unitsA.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
        if mouse in REVERSED:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unitsA.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
    
    temp_unitsB=[]
    for mouse in mice:
        if mouse in ORIGINAL:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unitsB.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
        if mouse in REVERSED:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_unitsB.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
                    
    plt.style.use('default')
    fig, ax=plt.subplots(1,1)                
    SEM1=sp.stats.sem(temp_unitsA, axis=0)
    SomHit, =ax.plot(np.arange(159), np.mean(temp_unitsA, axis=0), color='#A1007D')
    ax.fill_between(np.arange(159),np.mean(temp_unitsA, axis=0)+SEM1, np.mean(temp_unitsA, axis=0)-SEM1, color='#A1007D', alpha=0.1)
                      
    SEM2=sp.stats.sem(temp_unitsB, axis=0)
    VisHit, =ax.plot(np.arange(159), np.mean(temp_unitsB, axis=0), color='#00D88A')
    ax.fill_between(np.arange(159),np.mean(temp_unitsB, axis=0)+SEM1, np.mean(temp_unitsB, axis=0)-SEM1, color='#00D88A', alpha=0.1)
      
    #ax[i-1].grid(False)
    #ax[i-1].set_facecolor('white')
    ax.set_xticks([-0.5,39.5,79.5, 119.5, 159.5])
    #ax[i-1].tick_params(length=6, bottom=True)
    ax.set_xticklabels(['-1','0', '1', '2','3'])
    ax.set_yticks([0.5, 0.54, 0.58, 0.62])
    ax.set_yticklabels(['0.5','0.54','0.58','0.62'])
    ax.vlines(39.5, 0.48,0.62)
    ax.hlines(0.5, 1,158, linestyles='dotted')
    ax.set_ylim(0.48,0.62) #based on the range when showing the responses to HITS
    ax.set_xlim(-0.5,159.5)
    ax.set_xlabel('Seconds')
    ax.set_ylabel('MeanAUC')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend([SomHit, VisHit],['Right','Left'])
    fig.suptitle('Mean AUC - Stim Aligned - Hit vs CR\n All units, n=' + str(len(temp_unitsA)))
    
    
    #Scatter meanAUC (500msec after onset) - 
    
    
    temp_units3=[] 
    for mouse in mice:
        if mouse in ORIGINAL:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]
                        temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units3.append(0.5)
        if mouse in REVERSED:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]
                        temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units3.append(0.5)                     
    temp_units4=[]
    for mouse in mice:
        if mouse in ORIGINAL:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]
                        temp_units4.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units4.append(0.5)
        if mouse in REVERSED:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]
                        temp_units4.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units4.append(0.5)
                         
    plt.style.use('default')
    fig, ax=plt.subplots(1,1, figsize=(8,8))
    ax.fill([0.3,1,1],[0.3,1,0.3], '#A1007D', alpha=0.2)
    ax.fill([0.3,0.3,1],[0.3,1,1], '#00D88A', alpha=0.3)
    ax.set_xlim(0.3,1)
    ax.set_ylim(0.3,1)
    ax.scatter(temp_units3, temp_units4, c='k',s=30, linewidths=0)
    #MARK EXAMPLES
    
    ax.scatter(temp_units3[example1], temp_units4[example1],s=100, color='r', linewidths=0) #169
    ax.scatter(temp_units3[example2], temp_units4[example2],s=100, color='r', linewidths=0) #258
    ######
    ax.plot([0.3,1],[0.3,1], color='k')
    ax.vlines(0.5, 0.3,1, linestyles='dotted')
    ax.hlines(0.5,0.3,1, linestyles='dotted')
    ax.set_xlabel('Mean AUC Touch')
    ax.set_ylabel('Mean AUC Visual')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    s,p=sp.stats.wilcoxon(temp_units3, temp_units4)
    s,p1=sp.stats.wilcoxon([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))] , [y for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])

    fig.suptitle('Mean AUC ('+str(BinNumber*25)+'msec) - Stim Aligned - Hit vs CR \n All, n='+str(len(temp_units3))+ ' wilcoxon p='+str(p)[0:5] +'\n Sig: n='+
                 str(len([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])) +' wilcoxon p=' + str(p1)[0:5])
    
    #numers
    Right_mean=np.mean(temp_units3)
    Right_sem=np.std(temp_units3)/np.sqrt(len(temp_units3))
    Left_mean=np.mean(temp_units4)
    Left_sem=np.std(temp_units4)/np.sqrt(len(temp_units4))
    plt.plot(Right_mean, Left_mean, color='b', marker='o')
    plt.vlines(Right_mean, Left_mean-Left_sem, Left_mean+Left_sem) #np.log(Mean_CR-STD_CR)
    plt.hlines(Left_mean, Right_mean-Right_sem, Right_mean+Right_sem) #np.log(Mean_Hit-STD_Hit)
    return  Right_mean, Right_sem, Left_mean, Left_sem, p, p1


###############################################################################
# ITNB_scatter
###############################################################################

def ITNB_scatter(mice, master_log, example1, example2):
    total_temp_R=[]
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
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
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
    
    total_temp_L=[]
    temp_units2_totallicks=[]
    temp_mousename2=[]
    temp_KStest=[] 
    count=0
    temp_units2=[]           
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
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
    
    fig,ax=plt.subplots(1,1,figsize=(8,8))   
    #ax.scatter([x*2 for x in TouchBlock_Right],[x*2 for x in TouchBlock_Left], s=5, color='cornflowerblue')
    #ax.scatter([x*2 for x in VisualBlock_Right],[x*2 for x in VisualBlock_Left], s=5, color='orange') 
    ax.scatter([np.log(x*5) for x in peak_Right],[np.log(x*5) for x in peak_Left], s=50, color='k', alpha=0.5, linewidths=0)
    ax.plot([np.log(0.1),np.log(85)],[np.log(0.1),np.log(85)], color='grey')
    #MARK EXAMPLES
    ax.scatter(np.log(peak_Right[example1]*5), np.log(peak_Left[example1]*5),s=100, color='r', linewidths=0)
    ax.scatter(np.log(peak_Right[example2]*5), np.log(peak_Left[example2]*5),s=100, color='r', linewidths=0)
    ######
    ax.fill([np.log(0.1),np.log(85),np.log(85)],[np.log(0.1),np.log(0.1),np.log(85)], 'grey', alpha=0.3)
    #s,p1=sp.stats.wilcoxon(TouchBlock_Right, TouchBlock_Left)
    #s,p2=sp.stats.wilcoxon(VisualBlock_Right, VisualBlock_Left)
    s,p3=sp.stats.wilcoxon(peak_Right, peak_Left)
    ax.set_xlabel('Spontaneous Right Lick Responses (spikes/sec)', size=16)
    ax.set_ylabel('Spontaneous Left Lick Response (spikes/sec)', size=16)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.set_title('Spontaneous Licks - wilcoxon -\nTouchBlockPref p<'+str(p1)[-4:] + '\n VisualBlockPref p<'+str(p2)[-4:] + '\nRest p<'+str(p3)[:5])
    ax.set_title('Spontaneous Licks - wilcoxon - p<'+str(p3)[:5])
    plt.sca(ax)
    plt.xticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.yticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.xlim(np.log(0.1),np.log(85))
    plt.ylim(np.log(0.1),np.log(85))
    
    #numers
    peak_Right_mean=np.mean([x*5 for x in peak_Right])
    peak_Right_sem=np.std([x*5 for x in peak_Right])/np.sqrt(len(peak_Right))
    peak_Left_mean=np.mean([x*5 for x in peak_Left])
    peak_Left_sem=np.std([x*5 for x in peak_Left])/np.sqrt(len(peak_Left))
    plt.plot(np.log(peak_Right_mean), np.log(peak_Left_mean), color='b', marker='o')
    plt.vlines(np.log(peak_Right_mean), np.log(peak_Left_mean-peak_Left_sem), np.log(peak_Left_mean+peak_Left_sem)) #np.log(Mean_CR-STD_CR)
    plt.hlines(np.log(peak_Left_mean), np.log(peak_Right_mean-peak_Right_sem), np.log(peak_Right_mean+peak_Right_sem)) #np.log(Mean_Hit-STD_Hit)
    
    return  peak_Right_mean, peak_Right_sem, peak_Left_mean, peak_Left_sem, testR, testL,  temp_units1, temp_units2, temp_units1_totallicks, temp_units2_totallicks, temp_unit_names

###############################################################################
# example_ITBNandAUC
###############################################################################

def example_ITBNandAUC(mice, master_log, example, temp_units1, temp_units2, temp_units1_totallicks, temp_units2_totallicks, temp_unit_names):
    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
    REVERSED=['Claustrum18','Claustrum23','Claustrum25'] 
    
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
    plt.ylabel('Normalized Activity', size=16)
    plt.xlabel('Time (seconds)', size=16)
    plt.xlim(0,100)
    plt.yticks([0,2,4,6],['0','2','4', '6'], size=14)

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
    ax3.set_ylabel('Trials')
    
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
    plt.ylabel('Normalized Activity', size=16)
    plt.xlabel('Time (s)', size=16)
    plt.xlim(0,100)
    plt.yticks([0,2,4,6],['0','2','4', '6'], size=14)
    
    return

#########################################################################
# R2_correlation
#########################################################################
def R2_correlation(mice, master_log, testR, testL, example1, example2):
    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
    REVERSED=['Claustrum18','Claustrum23','Claustrum25'] 
    ITNB_meanDiff=[np.sum(testR[i][76:84])-np.sum(testL[i][76:84]) for i,x in enumerate(testR)]
    
    #HITs (lick aligned)
    temp_units1=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                norm=np.mean(np.mean(unit_log['-2to2sec_25msecbins_LickAligned'],0)[:39])
                temp_trial=[]
                if mouse in ORIGINAL:
                    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
                    if norm:
                        temp_units1.append(np.mean(temp_trial, axis=0))
                    else:
                        temp_units1.append(np.mean(temp_trial, axis=0))
                if mouse in REVERSED:
                    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
                    if norm:
                        temp_units1.append(np.mean(temp_trial, axis=0))
                    else:
                        temp_units1.append(np.mean(temp_trial, axis=0))
    temp_units2=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                norm=np.mean(np.mean(unit_log['-2to2sec_25msecbins_LickAligned'],0)[:39])
                temp_trial=[]
                if mouse in ORIGINAL:
                    for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
                    if norm:
                        temp_units2.append(np.mean(temp_trial, axis=0))
                    else:
                        temp_units2.append(np.mean(temp_trial, axis=0))
                if mouse in REVERSED:
                    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_trial.append(trial_log['-2to2sec_25msecbins_LickAligned']) 
                    if norm:
                        temp_units2.append(np.mean(temp_trial, axis=0))
                    else:
                        temp_units2.append(np.mean(temp_trial, axis=0))
           
    
    HIT_meanDiff=[np.sum(temp_units1[i][76:84])-np.sum(temp_units2[i][76:84]) for i,x in enumerate(temp_units1)] #reverse 1 and 2 depending on ORIGNAL or REVERSED
    
    
    #least-squares line
    X=np.asarray([x*5 for x in ITNB_meanDiff])
    y=np.asarray([x*5 for x in HIT_meanDiff])
    denominator=X.dot(X) - X.mean()*X.sum()
    m=(X.dot(y)-y.mean()*X.sum())/denominator
    b=(y.mean()*X.dot(X)-X.mean()*X.dot(y))/denominator
    y_pred=m*X+b
    
    res=y-y_pred
    tot=y-y.mean()
    R_squared=1 - res.dot(res)/tot.dot(tot)
    R=np.sqrt(R_squared)
    r,p=sp.stats.pearsonr(X,y)
    sp.stats.spearmanr(X,y)
    
    fig, ax=plt.subplots(1,1, figsize=(8,8))
    ax.set_yticks([-30, -15, 0, 15, 30])
    ax.set_yticklabels(['-30','-15','0','15','30'], size=14)
    ax.set_xticks([-30, -15, 0, 15, 30])
    ax.set_xticklabels(['-30','-15','0','15','30'], size=14)
    plt.sca(ax)
    linefit=plt.plot(X,y_pred, 'k')
    plt.legend(['R2:'+str(R_squared)[:4]+', p<'+str(p)])
    #points=plt.scatter([np.log(x) if x>0 else -(np.log(abs(x))) for x in X], [np.log(x) if x>0 else -(np.log(abs(x))) for x in y], c='k', alpha=0.2, linewidths=0, s=50)
    points=plt.scatter(X, y, c='k', alpha=0.2, linewidths=0, s=50)
    plt.xlim(-45,45)
    plt.ylim(-45,45)
    
    plt.xlabel('ITNB_meanDiff', size=16)
    plt.ylabel('HIT_meanDiff', size=16)
    plt.title('All \n Correlation betwwen ITNB-lick side preference \n and Hit-lick side preference')
    
    plt.scatter(X[example1], y[example1], s=20, color='r')
    plt.scatter(X[example2], y[example2], s=20, color='r')
    
    return

###############################################################################
# Loglikelihood_decoder
###############################################################################

def Loglikelihood_decoder(mice, master_log):
    score=[]
    Number_of_units=[]
    mouse_list=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            #First identify if day is usable: more than 3 units with significant response to ITNB
            usable_units=[]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit)
                    p1=1
                    p2=1
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    lick_num1=0
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        for lick in trial_log['NR_Right_licks_aligned_spikes_NonBouts']:
                            lick_num1+=1
                            temp_pre.append(np.sum((lick>-2) & (lick<-1.5)))
                            temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
    
                    K,p1=stats.ks_2samp(temp_pre,temp_post)
                
                    
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    lick_num2=0
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        for lick in trial_log['NR_Left_licks_aligned_spikes_NonBouts']:
                            lick_num2+=1
                            temp_pre.append(np.sum((lick>-2) & (lick<-1.5)))
                            temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                    
                    K,p2=stats.ks_2samp(temp_pre,temp_post)
                    
                    if (p1<0.05) | (p2<0.05):
                        if (lick_num1>30) & (lick_num2>30):
                            usable_units.append(1)
                        else:
                            usable_units.append(0)
                    else:
                        usable_units.append(0)
                        
    
            if np.sum(usable_units)>5:
                
                Number_of_units.append(usable_units)
                #First count the number of trials you will be using (SomHits and VisHits)
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], np.unique(day_log['cluster_name'])[0])]
                mouse_list.append(unit_log['mouse_name'].values[0])
                #Need to use the same amount of trial for each type, so find the smallest amount and use that
                Relevant_trials=np.min([len(unit_log[unit_log['Stim/Block/Response']=='VisHit']) , len(unit_log[unit_log['Stim/Block/Response']=='SomHit']) ] )
                #build matrices of the right size
                TuningMatrix=np.zeros((np.sum(usable_units),2))
                PopMatrix=np.zeros((  2*Relevant_trials  ,  np.sum(usable_units) ))
                TruthMatrix=np.zeros(( 2* Relevant_trials , 1))
                
                unit_counter=0
                usable_counter=0
                for unit in np.unique(day_log['cluster_name']):
                    if usable_units[usable_counter]:
                        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                        print(unit)
                        
                        p1=1
                        p2=1
                        temp_post=[] # to store dist post
                        lick_num=0
                        for trial in unit_log.index:
                            trial_log=unit_log.loc[trial,:]
                            for lick in trial_log['NR_Right_licks_aligned_spikes_NonBouts']:
                                lick_num+=1
                                temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                        R_tuningcurve=np.sum(temp_post)/lick_num
                        TuningMatrix[unit_counter,1]=R_tuningcurve
                    
                        
                        temp_post=[] # to store dist post
                        lick_num=0
                        for trial in unit_log.index:
                            trial_log=unit_log.loc[trial,:]
                            for lick in trial_log['NR_Left_licks_aligned_spikes_NonBouts']:
                                lick_num+=1
                                temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                        L_tuningcurve=np.sum(temp_post)/lick_num
                        TuningMatrix[unit_counter,0]=L_tuningcurve
                        
                        if mouse in ['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']: 
                            trial_counter=0    
                            for trial in unit_log.index:
                                if trial_counter<2*Relevant_trials:
                                     trial_log=unit_log.loc[trial,:]
                                     if trial_log['Stim/Block/Response']=='SomHit':
                                        data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                        PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                        TruthMatrix[trial_counter]=1
                                        trial_counter+=1
                                     if trial_log['Stim/Block/Response']=='VisHit':
                                        data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                        PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                        TruthMatrix[trial_counter]=0
                                        trial_counter+=1
                        elif mouse in ['Claustrum18', 'Claustrum23', 'Claustrum25']: 
                            trial_counter=0    
                            for trial in unit_log.index:
                                if trial_counter<2*Relevant_trials:
                                     trial_log=unit_log.loc[trial,:]
                                     if trial_log['Stim/Block/Response']=='SomHit':
                                        data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                        PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                        TruthMatrix[trial_counter]=0
                                        trial_counter+=1
                                     if trial_log['Stim/Block/Response']=='VisHit':
                                        data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                        PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                        TruthMatrix[trial_counter]=1
                                        trial_counter+=1
                        unit_counter+=1
                    usable_counter+=1
                #remove lick bias from tuning curves:
                Log_TuningMatrix=np.log(TuningMatrix)
                Log_TuningMatrix[:,0]= Log_TuningMatrix[:,0]-np.mean( Log_TuningMatrix[:,0])
                Log_TuningMatrix[:,1]= Log_TuningMatrix[:,1]-np.mean( Log_TuningMatrix[:,1])
                print(Log_TuningMatrix)
                Log_likelihood=np.matmul(PopMatrix, Log_TuningMatrix)  
                prediction=[]
                for x in Log_likelihood:
                    if x[1]>x[0]:
                        prediction.append(1)
                    elif x[1]<x[0]:
                        prediction.append(0)
                    else:
                        prediction.append(float('Nan'))
                Correct_predictions=[1  for i,x in enumerate(prediction) if x==TruthMatrix[i]]
                score.append((np.sum(Correct_predictions)/(2*Relevant_trials)))
                    
                    
    
    ORIGINAL_score=np.array(score)                
    
    
    score=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            #First identify if day is usable: more than 3 units with significant response to ITNB
            usable_units=[]
            for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    #norm=np.sum(np.mean(unit_log['-2to2sec_25msecbins_LickAligned'],0)[:39]) #baseline number of spikes per 1sec
                    print(unit)
                    p1=1
                    p2=1
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    lick_num1=0
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        for lick in trial_log['NR_Right_licks_aligned_spikes_NonBouts']:
                            lick_num1+=1
                            temp_pre.append(np.sum((lick>-2) & (lick<-1.5)))
                            temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                    
                    K,p1=stats.ks_2samp(temp_pre,temp_post)
                
                    
                    temp_pre=[] # to store dist pre
                    temp_post=[] # to store dist post
                    lick_num2=0
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        for lick in trial_log['NR_Left_licks_aligned_spikes_NonBouts']:
                            lick_num2+=1
                            temp_pre.append(np.sum((lick>-2) & (lick<-1.5)))
                            temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                    
                    K,p2=stats.ks_2samp(temp_pre,temp_post)
                    
                    if (p1<0.05) | (p2<0.05):
                        if (lick_num1>30) & (lick_num2>30):
                            usable_units.append(1)
                        else:
                            usable_units.append(0)
                    else:
                        usable_units.append(0)
                        
    
            if np.sum(usable_units)>5:
                #First count the number of trials you will be using (SomHits and VisHits)
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], np.unique(day_log['cluster_name'])[0])]
                #Need to use the same amount of trial for each type, so find the smallest amount and use that
                Relevant_trials=np.min([len(unit_log[unit_log['Stim/Block/Response']=='VisHit']) , len(unit_log[unit_log['Stim/Block/Response']=='SomHit']) ] )
                #build matrices of the right size
                TuningMatrix=np.zeros((np.sum(usable_units),2))
                PopMatrix=np.zeros((  2*Relevant_trials  ,  np.sum(usable_units) ))
                TruthMatrix=np.zeros((  2*Relevant_trials , 1))
                
                unit_counter=0
                usable_counter=0
                for unit in np.unique(day_log['cluster_name']):
                    if usable_units[usable_counter]:
                        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                        #norm=np.sum(np.mean(unit_log['-2to2sec_25msecbins_LickAligned'],0)[:39]) #baseline number of spikes per 1sec
                        print(unit)
                        #here need to build a random pseudo-session that has same amount of both type of trials
                        Som_list=np.random.choice(unit_log[unit_log['Stim/Block/Response']=='SomHit'].index.values,Relevant_trials, replace=False) #select the as many trials as there are relevant_trials
                        Vis_list=np.random.choice(unit_log[unit_log['Stim/Block/Response']=='VisHit'].index.values,Relevant_trials, replace=False) #select the as many trials as there are relevant_trials
                        All_list=np.hstack((Som_list, Vis_list)) #these represent a 'pseudo session' with equal trialtype numbers (the larger type is sampled randomly, the smaller has all trials)
                        A_list=np.random.choice(All_list,Relevant_trials, replace=False) #select the as many trials as there are relevant_trials
                        B_list=[x for x in All_list if x not in A_list]
                        temp_post=[] # to store dist post
                        lick_num=0
                        for trial in unit_log.index.values: #no need to limit the trilas you use for tuning
                            trial_log=unit_log.loc[trial,:]
                            for lick in trial_log['NR_Right_licks_aligned_spikes_NonBouts']:
                                lick_num+=1
                                temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                        Right_tuningcurve=np.sum(temp_post)/lick_num
                        TuningMatrix[unit_counter,1]=Right_tuningcurve
                    
                        temp_post=[] # to store dist post
                        lick_num=0
                        for trial in unit_log.index.values: #no need to limit the trilas you use for tuning
                            trial_log=unit_log.loc[trial,:]
                            for lick in trial_log['NR_Left_licks_aligned_spikes_NonBouts']:
                                lick_num+=1
                                temp_post.append(np.sum((lick>-0.25) & (lick<0.25)))#use these to test the responsiveness but alo to build the tuning curve (ie. the mean pre-motor response)
                        Left_tuningcurve=np.sum(temp_post)/lick_num
                        TuningMatrix[unit_counter,0]=Left_tuningcurve
                        
                        
                        trial_counter=0    
                        for trial in unit_log.index:
                             trial_log=unit_log.loc[trial,:]
                             if trial in A_list:
                                data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                TruthMatrix[trial_counter]=1
                                trial_counter+=1
                             if trial in B_list:
                                data=trial_log['Zscored_-2to2sec_25msecbins_LickAligned'][70:90]
                                PopMatrix[trial_counter,unit_counter]=np.sum(data)
                                TruthMatrix[trial_counter]=0
                                trial_counter+=1
                        unit_counter+=1
                    usable_counter+=1
                #remove lick bias from tuning curves:
                Log_TuningMatrix=np.log(TuningMatrix)
                Log_TuningMatrix[:,0]= Log_TuningMatrix[:,0]-np.mean( Log_TuningMatrix[:,0])
                Log_TuningMatrix[:,1]= Log_TuningMatrix[:,1]-np.mean( Log_TuningMatrix[:,1])
                print(Log_TuningMatrix)
                Log_likelihood=np.matmul(PopMatrix, Log_TuningMatrix)  
                #Chance=np.max(((len(TruthMatrix)-np.sum(TruthMatrix))/len(TruthMatrix), np.sum(TruthMatrix)/len(TruthMatrix)))
                prediction=[]
                for x in Log_likelihood:
                    if x[1]>x[0]:
                        prediction.append(1)
                    elif x[1]<x[0]:
                        prediction.append(0)
                    else:
                        prediction.append(float('Nan'))
                Correct_predictions=[1  for i,x in enumerate(prediction) if x==TruthMatrix[i]]
                score.append((np.sum(Correct_predictions)/(2*Relevant_trials)))
                if score==0:
                    print(mouse+day)
                
    CONTROL_score=np.array(score)
    
    
                
    ##if reversed
    #CONTROL_score=1-CONTROL_score
    #ORIGINAL_score=1-ORIGINAL_score
    
    
    fig,ax=plt.subplots(1,1,figsize=(3,9))
    plt.sca(ax)
    plt.scatter(np.zeros_like(ORIGINAL_score)+np.random.random(len(ORIGINAL_score)), ORIGINAL_score, alpha=0.3, c='k', linewidths=0, s=80)
    plt.scatter(0.5,np.mean(ORIGINAL_score), c='r', linewidths=0, s=80)
    plt.vlines(0.5, np.mean(ORIGINAL_score)-(np.std(ORIGINAL_score)/np.sqrt(len(ORIGINAL_score)))   , np.mean(ORIGINAL_score)+(np.std(ORIGINAL_score)/np.sqrt(len(ORIGINAL_score))) , color='r' )
    s,p1=stats.wilcoxon([x-0.5 for x in ORIGINAL_score])
    plt.text(0.5, 0.2, str(p1))
    
    
    plt.scatter(np.zeros_like(CONTROL_score)+np.random.random(len(CONTROL_score))+10, CONTROL_score, alpha=0.3, c='k', linewidths=0, s=80)
    plt.scatter(10.5,np.mean(CONTROL_score), c='r', linewidths=0, s=80)
    plt.vlines(10.5, np.mean(CONTROL_score)-(np.std(CONTROL_score)/np.sqrt(len(CONTROL_score)))   , np.mean(CONTROL_score)+(np.std(CONTROL_score)/np.sqrt(len(CONTROL_score))) , color='r' )
    s,p2=stats.wilcoxon([x-0.5 for x in CONTROL_score])
    plt.text(10.5, 0.2, str(p2))
    
    plt.hlines(0.5,-1,12,linestyles='dotted')
    
    plt.ylim(0,1)
    plt.xlim(-1,12)
    plt.xticks([0.5, 10.5], ['ORIGINAL', 'CONTROL'], size=16)
    plt.yticks([0,0.25, 0.5,0.75, 1], ['0','0.25',',0.5','0.75','1'], size=16)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    s,p0=stats.wilcoxon(ORIGINAL_score, CONTROL_score) 
    plt.title('log(Likelihood) decoder \nn= ' + str(len(ORIGINAL_score)) + ' sessions, Wilcoxon:'+str(p0))
    
    Score_mean=np.mean(ORIGINAL_score)
    Score_sem=np.std(ORIGINAL_score)/np.sqrt(len(ORIGINAL_score))
    Control_mean=np.mean(CONTROL_score)
    Control_sem=np.std(CONTROL_score)/np.sqrt(len(CONTROL_score))
    return Score_mean, Score_sem, Control_mean, Control_sem, p0, p1, p2, Number_of_units