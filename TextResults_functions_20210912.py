###############################################################################
# Results in text
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
    # - Responsiveness_S1()
    # - SomCR_vs_VisCR()
    # - Num_of_licks()
    # - Corr_lick_FR_AND_HitvsFApreLick()
    # - R2_optoOnly()
    # - FRpre_Claustrum_vs_S1()
    # - XCorr_VisualBlock()
    # - XCorr_TouchBlock()
    # - XCorr_CorrectTrial()
    # - XCorr_IncorrectTrial()
    # - XCorr_strength_byLickSide()
    
###############################################################################
# S1 responsiveness
###############################################################################
def Responsiveness_S1(trial_type, S1master_log, mice):
        
    unit_names=[]
    p_values=[]
    Direction=[]
    for mouse in mice:
        mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
        print(mouse)
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            print(day)
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                #print(unit_log['unit_name'].values[0])
                temp_PRE=[]
                temp_POST=[]
                for trial in unit_log[(unit_log['Stim/Block/Response']==trial_type) ].index:
                    trial_log=unit_log.loc[trial,:]
                    temp_PRE.append(np.sum((trial_log['StimALigned_spike_times']<0)&(trial_log['StimALigned_spike_times']>-0.1)))
                    temp_POST.append(np.sum((trial_log['StimALigned_spike_times']>0)&(trial_log['StimALigned_spike_times']<0.1)))
                try:
                    s,p=sp.stats.ks_2samp(temp_PRE, temp_POST)
                except:
                    continue
                Direction.append(np.mean(temp_POST)-np.mean(temp_PRE))
                unit_names.append(trial_log['unit_name'])
                p_values.append(p)
                

    Touch_responsive_units=[]
    for p,name,d in zip(p_values, unit_names,Direction):
        if p<0.05:
            Touch_responsive_units.append(name)

    
    Number_Resp= len(Touch_responsive_units)
    Number_Total=len(unit_names)
    
    return Number_Resp, Number_Total

###############################################################################
# How many neurons distinguish SomCR vs VisCR based on sensory response
###############################################################################
    
def SomCR_vs_VisCR(master_log, mice, BinNumber):

    temp_units3=[] 
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomCR_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomCR_vs_VisCR'].values[0]
                        if ((Onset_bin>38)&(Onset_bin<60)):
                            temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomCR_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                            plt.figure()
                            plt.plot(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomCR_vs_VisCR'].values[0])
                            plt.title(unit_log['unit_name'].values[0])
                            print(unit_log['unit_name'].values[0])
    
    Num=len(temp_units3) #34 neurons can distinguish SomCR_vs_VisCR
    return Num


###############################################################################
# How many neurons distinguish SomCR vs VisCR based on sensory response
###############################################################################
def Num_of_licks(mice, master_log):
    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
    REVERSED=['Claustrum18', 'Claustrum23', 'Claustrum25'] 
    plt.style.use('default')
    
    temp_units1=[] #for number of licks
    ILIs_HIT=[] #for inter lick inetrvals
    trialtypes=['SomHit','VisHit']
    for trialtype in trialtypes:
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            if ((mouse in ORIGINAL) & (trialtype=='SomHit')):
                LickSide='licks_right'
            elif ((mouse in ORIGINAL) & (trialtype=='VisHit')):
                LickSide='licks_left'
            elif ((mouse in REVERSED) & (trialtype=='SomHit')):
                LickSide='licks_left'
            elif ((mouse in REVERSED) & (trialtype=='VisHit')):
                LickSide='licks_right'
            for day in np.unique(mouse_log['date']):
                print(mouse + day[0])
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    for trial in unit_log[unit_log['Stim/Block/Response']==trialtype].index:
                        trial_log=unit_log.loc[trial,:]
                        temp_units1.append(sum(trial_log[LickSide][0][0]<(trial_log['stim_onset'][0][0][0][0]+1))[0]) #this outputs the number of licks within 1 second of the first lick post stim
                        Hit_licks=trial_log[LickSide][0][0][trial_log[LickSide][0][0]<(trial_log['stim_onset'][0][0][0][0]+1)]
                        if len(Hit_licks)>=2:
                            ILIs_HIT.append(np.diff( Hit_licks))
    
     
    temp_units2=[]
    ILIs_FA=[] #for inter lick inetrvals
    trialtypes=['SomFA','VisFA']
    for trialtype in trialtypes:
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(mouse + day[0])
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    for trial in unit_log[unit_log['Stim/Block/Response']==trialtype].index:
                        trial_log=unit_log.loc[trial,:]
                        if isinstance((trial_log['FirstLick']), float):
                            print(trial_log['unit_name'])
                            continue
                        if len(trial_log['licks_right'][0][0])>0:
                            RightLicks=trial_log['licks_right'][0][0]
                            RightLicks=[x[0] for x in RightLicks]
                        else:
                            RightLicks=[]
                        if len(trial_log['licks_left'][0][0])>0:
                            LeftLicks=trial_log['licks_left'][0][0]
                            LeftLicks=[x[0] for x in LeftLicks]
                        else:
                            LeftLicks=[] 
                        Licks=RightLicks+LeftLicks
                        temp_units2.append(sum(Licks<(trial_log['stim_onset'][0][0][0][0]+1))) #this outputs the number of licks within 1 second of the first lick post stim
                        FA_licks=trial_log[LickSide][0][0][trial_log[LickSide][0][0]<(trial_log['stim_onset'][0][0][0][0]+1)]
                        if len(FA_licks)>=2:
                            ILIs_FA.append(np.diff( FA_licks))
                            
    
    fig,ax=plt.subplots(1,1,figsize=(5,5))
    Cum=stats.cumfreq(temp_units1, numbins=50, defaultreallimits=(0,50))#, defaultreallimits=(-0.5,1))
    x=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    x = [np.log(a) if a!=0 else 0 for a in x ]
    HITplot=plt.plot(x,Cum.cumcount/np.size(temp_units1), color='k')
    Cum=stats.cumfreq(temp_units2, numbins=50, defaultreallimits=(0,50))#, defaultreallimits=(-0.5,1))
    x=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    x = [np.log(a) if a!=0 else 0 for a in x ]
    FAplot=plt.plot(x,Cum.cumcount/np.size(temp_units2), color='k', linestyle='dotted')
    
    x_values=[np.log(x) for x in [5,10,20,40]]
    plt.xticks(x_values, ['5','10','20','40'])
    plt.xlabel('Number of licks')
    plt.ylabel('CDF')
    s,p=stats.mannwhitneyu(temp_units1, temp_units2)
    plt.legend(['HITs','FAs'])
    plt.title('Mann-WhitneyU p=' +str(p) +'\n'+ str(len(temp_units1))+' HIT trials, '+str(len(temp_units2))+' FA trials')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    meanHit=np.mean(temp_units1)
    semHit=np.std(temp_units1)/np.sqrt(len(temp_units1))
    nHit=len(temp_units1)
    meanFA=np.mean(temp_units2)
    semFA=np.std(temp_units2)/np.sqrt(len(temp_units2))
    nFA=len(temp_units2)
    #meanHit:7.20 ,semHit:0.0233 , meanFA:6.092 , semFA:0.0389 , nHit=71454 trials, nFA=24521 trials, p=3.0e-244    
    
    return meanHit, semHit, meanFA, semFA, nHit, nFA, nFA, p


###############################################################################
# Correlations between Number of licks and FR
###############################################################################
def Corr_lick_FR_AND_HitvsFApreLick(master_log, mice):    
    Direction=[]
    Touch_responsive_units=[]
    p_values_suround=[]
    p_values_pre=[]
    p_values_post=[]
    Hit_means_suround=[]
    FA_means_suround=[]
    Hit_means_pre=[]
    FA_means_pre=[]
    Hit_means_post=[]
    FA_means_post=[]
    Hit_licks=[]
    Hit_spikes_pre=[]
    Hit_spikes_post=[]
    Hit_spikes_suround=[]
    Bout_licks=[]
    Bout_length=[]
    Bout_rate=[]
    Hit_RT=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        print(mouse)
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            print(day)
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                #print(unit_log['unit_name'].values[0])
                temp_PRE=[]
                temp_POST=[]
                for trial in unit_log[(unit_log['Stim/Block/Response']=='SomHit') ].index:
                    trial_log=unit_log.loc[trial,:]
                    temp_PRE.append(np.sum((trial_log['StimALigned_spike_times']<0)&(trial_log['StimALigned_spike_times']>-1)))
                    temp_POST.append(np.sum((trial_log['StimALigned_spike_times']>0)&(trial_log['StimALigned_spike_times']<1)))
                s,p=sp.stats.ks_2samp(temp_PRE, temp_POST)
                if p<0.05:
                    Direction.append(np.mean(temp_POST)-np.mean(temp_PRE))
                    Touch_responsive_units.append(trial_log['unit_name'])
                    
                    temp_units1_spikes=[]
                    temp_units3_spikes=[]
                    temp_units5_spikes=[]
                    temp_units1_licks=[]
                    temp_units1_spikes_pre=[]
                    temp_units1_spikes_post=[]
                    temp_Bout_licks=[]
                    temp_Bout_length=[]
                    temp_Bout_rate=[]
                    RT=[]
                    for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                        trial_log=unit_log.loc[trial,:]
                        
                            
                        if isinstance((trial_log['FirstLick']), float):
                            print(trial_log['unit_name'])
                            continue
                        if len(trial_log['licks_right'][0][0])>0:
                            RightLicks=trial_log['licks_right'][0][0]
                            RightLicks=[x[0] for x in RightLicks]
                        else:
                            RightLicks=[]
                        if len(trial_log['licks_left'][0][0])>0:
                            LeftLicks=trial_log['licks_left'][0][0]
                            LeftLicks=[x[0] for x in LeftLicks]
                        else:
                            LeftLicks=[] 
                        Licks=RightLicks+LeftLicks
        
                        temp_units1_licks.append(sum((Licks<(trial_log['stim_onset'][0][0][0][0]+3)) & (Licks>(trial_log['stim_onset'][0][0][0][0])))) #this outputs the number of licks within 1 second of the first lick post stim
                        RT.append(trial_log['FirstLick'][0][0])
                        
                        lick_mask=(Licks<(trial_log['stim_onset'][0][0][0][0]+3)) & (Licks>(trial_log['stim_onset'][0][0][0][0]))
                        Relevant_licks=[x for x,m in zip(Licks, lick_mask) if m]
                        last_bout_lick_index=1 #set as default
                        if len(Relevant_licks)==1:
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>-1)))[0])
                                temp_units3_spikes.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-trial_log['FirstLick'][0][0])))[0] / trial_log['FirstLick'][0][0])
                                temp_units5_spikes.append(1)
                            else:
                                temp_units1_spikes.append(0)
                                temp_units3_spikes.append(0)
                                temp_units5_spikes.append(0)
                                
                            temp_Bout_licks.append(1)
                            temp_Bout_length.append(1)
                            temp_Bout_rate.append(1)          
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes_pre.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-1)))[0])
                            else:
                                temp_units1_spikes_pre.append(0)
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes_post.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>0)))[0])
                            else:
                                temp_units1_spikes_post.append(0)
                                
                        else:
                            ILI=np.diff(Relevant_licks)
                            ILI_mask=ILI<0.5
                            i=0
                            while (ILI_mask[i]) :
                                last_bout_lick_index=i+1
                                i+=1
                                if i==len(ILI_mask):
                                    break
     
                            temp_Bout_licks.append(last_bout_lick_index+1)
                            temp_Bout_length.append(sum(ILI[:last_bout_lick_index]))
                            if sum(ILI[:last_bout_lick_index])==0:
                                temp_Bout_rate.append(0)
                            else:
                                temp_Bout_rate.append(last_bout_lick_index+1/sum(ILI[:last_bout_lick_index]))
                                
                                    
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes_pre.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-1)))[0])
                            else:
                                temp_units1_spikes_pre.append(0)
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes_post.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>0)))[0])
                            else:
                                temp_units1_spikes_post.append(0)
                                
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units1_spikes.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>-1)))[0])
                                temp_units3_spikes.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-trial_log['FirstLick'][0][0])))[0] / trial_log['FirstLick'][0][0])
                                temp_units5_spikes.append(sum(((trial_log['LickALigned_spike_times']<sum(ILI[:last_bout_lick_index])) & (trial_log['LickALigned_spike_times']>0)))[0] / sum(ILI[:last_bout_lick_index]))
                            else:
                                temp_units1_spikes.append(0)
                                temp_units3_spikes.append(0)
                                temp_units5_spikes.append(0)
                            
                    Hit_RT.append(RT)
                    Hit_licks.append(temp_units1_licks)
                    Hit_spikes_pre.append(temp_units1_spikes_pre)
                    Hit_spikes_post.append(temp_units1_spikes_post)
                    Hit_spikes_suround.append(temp_units1_spikes)
                    Bout_licks.append( temp_Bout_licks)
                    Bout_length.append( temp_Bout_length)
                    Bout_rate.append( temp_Bout_rate)
                    
    #                r,p=stats.pearsonr(temp_units1_licks, temp_units1_spikes_post)
    #                if p>0.05:
    #                    fig,ax=plt.subplots(1,1,figsize=(5,5))
    #                    plt.scatter(temp_units1_licks, temp_units1_spikes_post)
    #                    plt.title(trial_log['unit_name'] + ' r='+ str(r))
                    
                    temp_units2_spikes=[]
                    temp_units4_spikes=[]
                    temp_units6_spikes=[]
                    for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
                        trial_log=unit_log.loc[trial,:]
                        
                        if isinstance((trial_log['FirstLick']), float):
                            print(trial_log['unit_name'])
                            continue
                        if len(trial_log['licks_right'][0][0])>0:
                            RightLicks=trial_log['licks_right'][0][0]
                            RightLicks=[x[0] for x in RightLicks]
                        else:
                            RightLicks=[]
                        if len(trial_log['licks_left'][0][0])>0:
                            LeftLicks=trial_log['licks_left'][0][0]
                            LeftLicks=[x[0] for x in LeftLicks]
                        else:
                            LeftLicks=[] 
                        Licks=RightLicks+LeftLicks
                            
                        lick_mask=(Licks<(trial_log['stim_onset'][0][0][0][0]+3)) & (Licks>(trial_log['stim_onset'][0][0][0][0]))
                        Relevant_licks=[x for x,m in zip(Licks, lick_mask) if m]
                        last_bout_lick_index=1 #set as default
                        if len(Relevant_licks)==1:
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units2_spikes.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>-1)))[0])
                                temp_units4_spikes.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-trial_log['FirstLick'][0][0])))[0] / trial_log['FirstLick'][0][0])
                                temp_units6_spikes.append(1)
                            else:
                                temp_units2_spikes.append(0)
                                temp_units4_spikes.append(0)
                                temp_units6_spikes.append(0)
                        else:
                            ILI=np.diff(Relevant_licks)
                            ILI_mask=ILI<0.5
                            i=0
                            while (ILI_mask[i]) :
                                last_bout_lick_index=i+1
                                i+=1
                                if i==len(ILI_mask):
                                    break
                            if len(trial_log['LickALigned_spike_times'])>0:
                                temp_units2_spikes.append(sum(((trial_log['LickALigned_spike_times']<1) & (trial_log['LickALigned_spike_times']>-1)))[0])
                                temp_units4_spikes.append(sum(((trial_log['LickALigned_spike_times']<0) & (trial_log['LickALigned_spike_times']>-trial_log['FirstLick'][0][0])))[0] / trial_log['FirstLick'][0][0])
                                temp_units6_spikes.append(sum(((trial_log['LickALigned_spike_times']<sum(ILI[:last_bout_lick_index])) & (trial_log['LickALigned_spike_times']>0)))[0] / sum(ILI[:last_bout_lick_index]))
                            else:
                                temp_units2_spikes.append(0)
                                temp_units4_spikes.append(0)
                                temp_units6_spikes.append(0)
                    s,p=stats.mannwhitneyu(temp_units1_spikes, temp_units2_spikes)
                    p_values_suround.append(p)
                    Hit_means_suround.append(np.mean(temp_units1_spikes))
                    FA_means_suround.append(np.mean(temp_units2_spikes))
                    
                    s,p=stats.mannwhitneyu(temp_units3_spikes, temp_units4_spikes)
                    p_values_pre.append(p)
                    Hit_means_pre.append(np.mean(temp_units3_spikes))
                    FA_means_pre.append(np.mean(temp_units4_spikes))
                    
                    if sum(temp_units5_spikes)>0:
                        s,p=stats.mannwhitneyu(temp_units5_spikes, temp_units6_spikes)
                        p_values_post.append(p)
                        Hit_means_post.append(np.mean(temp_units5_spikes))
                        FA_means_post.append(np.mean(temp_units6_spikes))
    
                
    spikes=Hit_spikes_suround
    UP=[]
    DOWN=[]
    UP_Corr_pos=0
    UP_Corr_neg=0
    DOWN_Corr_pos=0
    DOWN_Corr_neg=0
    for a,b,c in zip(Bout_licks, spikes, Direction):
        r,p=stats.pearsonr(a,b)
        if p<0.05:
            if c<0:
                DOWN.append(r)
                if r>0:
                    DOWN_Corr_pos+=1
                if r<0:
                    DOWN_Corr_neg+=1
            if c>0:
                UP.append(r)
                if r>0:
                    UP_Corr_pos+=1
                if r<0:
                    UP_Corr_neg+=1
            
    plt.figure()        
    plt.hist(UP)
    plt.hist(DOWN)
    plt.title(str(len(UP)) +' activated neurons ('+str(UP_Corr_pos)+' pos, ' + str(UP_Corr_neg) + ' neg), \n'
              +str(len(DOWN))+ ' inhibited neurons('+str(DOWN_Corr_pos)+' pos, ' + str(DOWN_Corr_neg) + 
              ') are correlated\nwith number of licks in response')
    
    
    #  75(60 activated (25pos_corr, 35 neg_corr) +15 inhibited (7pos_cor, 8neg_cor)) out of 347 responsive (21.6%)

    #Numbers to show that activity is stronger in Hits even before the reward comes on
    s,p_inh=stats.wilcoxon([x for x,c in zip(Hit_means_pre,Direction) if c<0], [x for x,c in zip(FA_means_pre,Direction) if c<0])
    s,p_exc=stats.wilcoxon([x for x,c in zip(Hit_means_pre,Direction) if c>0], [x for x,c in zip(FA_means_pre,Direction) if c>0])
    Hitmean_exc=np.mean([x for x,c in zip(Hit_means_pre,Direction) if c>0] )
    Hitmsem_exc=np.std([x for x,c in zip(Hit_means_pre,Direction) if c>0] )/np.sqrt(len([x for x,c in zip(Hit_means_pre,Direction) if c>0] ))
    FAmean_exc=np.mean([x for x,c in zip(FA_means_pre,Direction) if c>0] )
    FAmsem_exc=np.std([x for x,c in zip(FA_means_pre,Direction) if c>0] )/np.sqrt(len([x for x,c in zip(FA_means_pre,Direction) if c>0] ))
    
    Hitmean_inh=np.mean([x for x,c in zip(Hit_means_pre,Direction) if c<0] )
    Hitmsem_inh=np.std([x for x,c in zip(Hit_means_pre,Direction) if c<0] )/np.sqrt(len([x for x,c in zip(Hit_means_pre,Direction) if c<0] ))
    FAmean_inh=np.mean([x for x,c in zip(FA_means_pre,Direction) if c<0] )
    FAmsem_inh=np.std([x for x,c in zip(FA_means_pre,Direction) if c<0] )/np.sqrt(len([x for x,c in zip(FA_means_pre,Direction) if c<0] ))
    
    #For inhibited neurons: p=4.576e-4, meanHit=4.645 ,semHit=0.527 , meanFA=5.387 ,semFA=0.585 , n=87
    #For activated neurons: p=2.550e-15, meanHit=15.50 ,semHit=0.926 , meanFA=13.328 ,semFA=0.8246 , n=260
    
    return Hitmean_exc, Hitmsem_exc,FAmean_exc,FAmsem_exc,Hitmean_inh, Hitmsem_inh, FAmean_inh, FAmsem_inh



########################################################################
# plot optoneurons on the R2 plot
########################################################################
def R2_optoOnly(mice, master_log, testR, testL):
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
    
    
    
    #Get the indices of opto in the lists tat are used for plotting
    OptoInd=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                if unit_log['Category'].values[0]=='OptoTag':
                    OptoInd.append(1)
                else:
                    OptoInd.append(0)
                    
                
    # only Opto on top
    #least-squares line
    X=np.asarray([x*5 for x,y in zip(ITNB_meanDiff,OptoInd) if y==1])
    y=np.asarray([x*5 for x,y in zip(HIT_meanDiff,OptoInd) if y==1])
    denominator=X.dot(X) - X.mean()*X.sum()
    m=(X.dot(y)-y.mean()*X.sum())/denominator
    b=(y.mean()*X.dot(X)-X.mean()*X.dot(y))/denominator
    y_pred=m*X+b
    
#    res=y-y_pred
#    tot=y-y.mean()
#    R_squared=1 - res.dot(res)/tot.dot(tot)
#    R=np.sqrt(R_squared)
#    r,p=sp.stats.pearsonr(X,y)
#    sp.stats.spearmanr(X,y)
#    
#    fig, ax=plt.subplots(1,1, figsize=(8,8))
#    ax.set_yticks([-30, -15, 0, 15, 30])
#    ax.set_yticklabels(['-30','-15','0','15','30'], size=14)
#    ax.set_xticks([-30, -15, 0, 15, 30])
#    ax.set_xticklabels(['-30','-15','0','15','30'], size=14)
#    plt.sca(ax)
#    linefit=plt.plot(X,y_pred, 'k')
#    plt.legend(['R2:'+str(R_squared)[:4]+', p<'+str(p)])
#    #points=plt.scatter([np.log(x) if x>0 else -(np.log(abs(x))) for x in X], [np.log(x) if x>0 else -(np.log(abs(x))) for x in y], c='k', alpha=0.2, linewidths=0, s=50)
#    points=plt.scatter(X, y, c='b', alpha=0.8, linewidths=0, s=30)
#    plt.xlim(-45,45)
#    plt.ylim(-45,45)
    return r,p


########################################################################
# plot optoneurons on the R2 plot
########################################################################
def FRpre_Claustrum_vs_S1(XCorr_df, S1XCorr_df):
    FRpre=[]
    for pair in XCorr_df.index: 
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if ((Centered_SigBins>0) & (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name'])):
            print(count)
    #        if ( (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name'])):
            Corr_LickPrefs1.append(XCorr_df.loc[pair,'Neuron1_LickPref'])
            Corr_LickPrefs2.append(XCorr_df.loc[pair,'Neuron2_LickPref'])
            Corr_strength.append( np.max(Subtracted_XCorr[488:512]) - Baseline )  
            FRpre.append(XCorr_df.loc[pair,'Geometric mean FR'])
    
    S1FRpre=[]     
    for pair in S1XCorr_df.index:
        Subtracted_XCorr=S1XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-S1XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        if ((Centered_SigBins>0) & (S1XCorr_df.loc[pair,'Neuron1_name'] != S1XCorr_df.loc[pair,'Neuron2_name'])):
    #        if ( (S1XCorr_df.loc[pair,'Neuron1_name'] != S1XCorr_df.loc[pair,'Neuron2_name'])):
            S1Corr_LickPrefs1.append(S1XCorr_df.loc[pair,'Neuron1_LickPref'])
            S1Corr_LickPrefs2.append(S1XCorr_df.loc[pair,'Neuron2_LickPref'])
            S1Corr_strength.append( np.max(Subtracted_XCorr[488:512]) - Baseline )       
            S1FRpre.append(S1XCorr_df.loc[pair,'Geometric mean FR'])
    plt.figure()
    sns.kdeplot(FRpre, color='b', alpha=0.5)
    sns.kdeplot(S1FRpre,  color='r', alpha=0.5)
    s,p=stats.ttest_ind(FRpre, S1FRpre)
    Mean=np.mean(FRpre)
    SEM=np.std(FRpre)/np.sqrt(len(FRpre))
    n=len(FRpre)
    S1Mean=np.mean(S1FRpre)
    S1SEM=np.std(S1FRpre)/np.sqrt(len(S1FRpre))
    nS1=len(S1FRpre)
    plt.title('FR pre, p='+str(p) + '\nClaustrum - mean:'+str(Mean)+', sem:'+str(SEM)+'\nS1 - mean:'+str(S1Mean)+', sem:'+str(S1SEM))
    #Results: meanFR:8.75 , semFR:0.674 , MeanFRS1:12.00 , semFRS1:0.639 , n=75 , nS1=129, p=0.00118
    return Mean, SEM, S1Mean, S1SEM, n, nS1, p

########################################################################
# XCorr proportions across various groups
########################################################################
def XCorr_VisualBlock(XCorr_df):
    total_count=0
    count=0
    #Sig_index=[]
    for pair in XCorr_df.index:
        
        if XCorr_df.loc[pair,'Neuron1_name'][-20:-18] in ['18','23','25']:
            Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_RIGHT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_RIGHT1sec']
        elif XCorr_df.loc[pair,'Neuron1_name'][-20:-18] in ['l4','l5','l6','31','32','37']:
            Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_LEFT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_LEFT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']):
            total_count+=1
            if (Centered_SigBins>0) :
                count+=1
    
    return count, total_count

def XCorr_TouchBlock(XCorr_df):
    total_count=0
    count=0
    #Sig_index=[]
    for pair in XCorr_df.index:
        
        if XCorr_df.loc[pair,'Neuron1_name'][-20:-18] in ['18','23','25']:
            Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_LEFT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_LEFT1sec']
        elif XCorr_df.loc[pair,'Neuron1_name'][-20:-18] in ['l4','l5','l6','31','32','37']:
            Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_RIGHT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_RIGHT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']):
            total_count+=1
            if (Centered_SigBins>0) :
                count+=1
    
    return count, total_count

def XCorr_CorrectTrial(XCorr_df):
    total_count=0
    count=0
    #Sig_index=[]
    for pair in XCorr_df.index:
        
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_CORRECT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_CORRECT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']):
            total_count+=1
            if (Centered_SigBins>0) :
                count+=1
    return count, total_count

def XCorr_IncorrectTrial(XCorr_df):
    total_count=0
    count=0
    #Sig_index=[]
    for pair in XCorr_df.index:
        
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_INCORRECT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_INCORRECT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']):
            total_count+=1
            if (Centered_SigBins>0) :
                count+=1
    return count, total_count

########################################################################
# Is strength of Right lick neurons that are corr the same during ipsi block or contra blocks
########################################################################
def XCorr_strength_byLickSide(XCorr_df):
    #get list of corr and list of Right_pref

    count=0
    Sig_index=[]
    for pair in XCorr_df.index:
        
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])        
        if ((Centered_SigBins>0) & (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name'])):
            Sig_index.append(pair)
        count+=1
    
    # now repeat Corr calculation but using RIGHT or LEFT, ad limiting to those neurons
    
    
    RIGHT_strength=[]
    for pair in Sig_index:
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_RIGHT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_RIGHT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if ((XCorr_df.loc[pair,'Neuron1_LickPref']>0.1) & (XCorr_df.loc[pair,'Neuron2_LickPref']>0.1)  & (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']) ):
            RIGHT_strength.append(np.max(Subtracted_XCorr[488:512]) - Baseline)
    
    
    LEFT_strength=[]
    for pair in Sig_index:
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_LEFT1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_LEFT1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        BaselineSD=np.std( np.hstack((Subtracted_XCorr[:100] , Subtracted_XCorr[-100:])))
        Top_cutoff=Baseline+5*BaselineSD
        Bottom_cutoff=Baseline-5*BaselineSD
        Sig_bins=[1 if ((x>Top_cutoff) | (x<Bottom_cutoff)) else 0 for x in  Subtracted_XCorr]
        Centered_SigBins=np.sum(Sig_bins[488:512])
        
        if ((XCorr_df.loc[pair,'Neuron1_LickPref']>0.1) & (XCorr_df.loc[pair,'Neuron2_LickPref']>0.1)  & (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name']) ):
            LEFT_strength.append(np.max(Subtracted_XCorr[488:512]) - Baseline)
            
    s,p=stats.wilcoxon( LEFT_strength,  RIGHT_strength)
    meanRIGHT=np.mean( RIGHT_strength)
    semRIGHT=np.std(RIGHT_strength)/np.sqrt(len(RIGHT_strength))
    meanLEFT=np.mean( LEFT_strength)
    semLEFT=np.std(LEFT_strength)/np.sqrt(len(LEFT_strength))
    return meanRIGHT, semRIGHT, meanLEFT, semLEFT, p


########################################################################
# ILI rate DREADD
########################################################################
def ILIrate_DREADD(mice, master_DREADD):
    
    def Lick_analysis(day_log, trial_range, trial_type='All'):
      RT=[]
      temp_rewardedlicks=[]
      temp_nonrewardedlicks=[]
      temp_trial_num=[]
      temp_rewardedlicks_trial_length=[]
      temp_nonrewardedlicks_trial_length=[]
      
      
        
      for trial in day_log.index[trial_range[0]:trial_range[1]]:
            trial_log=day_log.loc[trial,:]
            if trial_type!='All':
                  if trial_log['Stim/Block/Response']!=trial_type:
                        continue
          
            # Fisrt get data related to rewarded licks
            if len(trial_log['rewarded_licks'])==1:
                  temp_rewardedlicks.append(trial_log['rewarded_licks'][0])    
                  temp_rewardedlicks_trial_length.append(trial_log['Trial_Start_Stop'][0][1] - trial_log['Trial_Start_Stop'][0][0])
            else:     
                  temp_rewardedlicks.append([x[0] for x in trial_log['rewarded_licks']])
                  temp_rewardedlicks_trial_length.append(trial_log['Trial_Start_Stop'][0][1] - trial_log['Trial_Start_Stop'][0][0])

            #Get data related to nonrewarded licks
            if len(trial_log['non_rewarded_licks'])>0:
                  temp_nonrewardedlicks.append(trial_log['non_rewarded_licks'])
                  temp_nonrewardedlicks_trial_length.append(trial_log['Trial_Start_Stop'][0][1] - trial_log['Trial_Start_Stop'][0][0])
                  
                  #get RT     
#            if trial_log['FirstLick']>1:
#                  continue
            if isinstance((trial_log['FirstLick']), float):
                  continue
            else:
                  RT.append(trial_log['FirstLick'][0][0])
            
            temp_trial_num.append(trial_log['trial_num'][0][0][0][0])
                  
      return RT, temp_rewardedlicks, temp_nonrewardedlicks, temp_trial_num, temp_rewardedlicks_trial_length, temp_nonrewardedlicks_trial_length



    ## Run with rate of ITNB 
    trial_type='All'
    trial_range=[0,-1]
    Summary_data=[]
    for virus in ['DREADD','CONTROL']:
          ALL_CONTROL_SALINE=[]
          ALL_CONTROL_AGONIST21=[]
          ALL_MICE=[]
          for mouse in mice:
                SALINE_data=[]
                AGONIST21_data=[]
                mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]
                if mouse_log['virus'].values[0]==virus:
                        print(mouse)
                        for day in np.unique(mouse_log['date'])[:8]:
                            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                            print(day)
                            RT, temp_rewardedlicks, temp_nonrewardedlicks, temp_trial_num, temp_rewardedlicks_trial_length, temp_nonrewardedlicks_trial_length = Lick_analysis(day_log, trial_range, trial_type)
                            temp_nonrewardedlicks=[len(x)/y for x,y in  zip(temp_nonrewardedlicks, temp_nonrewardedlicks_trial_length)]
                            if day_log['Injection'].values[0]=='SALINE':
                                  SALINE_data.append([np.mean(temp_nonrewardedlicks)])
                            if day_log['Injection'].values[0]=='AGONIST21':
                                  AGONIST21_data.append([np.mean(temp_nonrewardedlicks)])
                        ALL_CONTROL_SALINE.append(SALINE_data)
                        ALL_CONTROL_AGONIST21.append(AGONIST21_data)
                        ALL_MICE.append(mouse)
          
          #test difference in fraction of a trial type
          All_Saline=[]
          All_Agonist21=[]
          All_mouse=[]
          XX=[]
          YY=[]
          ZZ=[]
          for mouse_saline,mouse_agonist,mouse in  zip(ALL_CONTROL_SALINE, ALL_CONTROL_AGONIST21, ALL_MICE):
                SALINE_data=[]
                AGONIST21_data=[]
                for a,b in zip(mouse_saline,mouse_agonist):
                      All_Saline.append(a[0]) 
                      All_Agonist21.append(b[0] )
                      SALINE_data.append(a[0]) 
                      AGONIST21_data.append(b[0]) 
                    
                      All_mouse.append(mouse)
                s,p=stats.mannwhitneyu(SALINE_data,AGONIST21_data)
                print('_pairedSessions_'+str(p))
                
                XX.append(np.mean(SALINE_data))
                YY.append(np.mean(AGONIST21_data))
                ZZ.append(mouse)
          s,p=stats.wilcoxon(XX, YY)
          print('pairedMice_'+str(p))
                      
          # Summary_data.append([All_Saline, All_Agonist21, All_mouse])
          Summary_data.append([XX,YY, ZZ])
          
    dreadd_plot=Summary_data[0]
    control_plot=Summary_data[1]
    
    colors=['b','r','g','k','orange','purple']
    fig,ax=plt.subplots(1,1,figsize=(5,15))
    for i,mouse in enumerate(np.unique(dreadd_plot[2])):
        data1=[x for x,y in zip(dreadd_plot[0], dreadd_plot[2]) if ( (y==mouse) & (x!=0) )]
        # plt.scatter(np.zeros_like(data1), data1, c='k', alpha=0.2, linewidths=0)
        data2=[x for x,y in zip(dreadd_plot[1], dreadd_plot[2]) if ( (y==mouse) & (x!=0) )]
        # plt.scatter(np.ones_like(data2), data2, c='k', alpha=0.2, linewidths=0)
        plt.plot([0,1], [np.mean(data1), np.mean(data2)], color='k', alpha=0.2)
    MEAN_DREADD_SALINE_acrossMice=np.mean(dreadd_plot[0])
    MEAN_DREADD_AGONIST_acrossMice=np.mean(dreadd_plot[1])
    SEM_DREADD_SALINE_acrossMice=np.std(dreadd_plot[0])/np.sqrt(len(dreadd_plot[0]))
    SEM_DREADD_AGONIST_acrossMice=np.std(dreadd_plot[1])/np.sqrt(len(dreadd_plot[1]))
    
    plt.plot([0,1], [MEAN_DREADD_SALINE_acrossMice, MEAN_DREADD_AGONIST_acrossMice], color='red', linewidth=4)
    plt.vlines(0,MEAN_DREADD_SALINE_acrossMice-SEM_DREADD_SALINE_acrossMice, MEAN_DREADD_SALINE_acrossMice+SEM_DREADD_SALINE_acrossMice, color='red', linewidths=4)
    plt.vlines(1,MEAN_DREADD_AGONIST_acrossMice-SEM_DREADD_AGONIST_acrossMice, MEAN_DREADD_AGONIST_acrossMice+SEM_DREADD_AGONIST_acrossMice, color='red', linewidths=4)
    
    
    for i,mouse in enumerate(np.unique(control_plot[2])):
        data1=[x for x,y in zip(control_plot[0], control_plot[2]) if ( (y==mouse) & (x!=0) )]
        # plt.scatter(np.zeros_like(data1)+2, data1, c='k', alpha=0.2, linewidths=0)
        data2=[x for x,y in zip(control_plot[1], control_plot[2]) if ( (y==mouse) & (x!=0) )]
        # plt.scatter(np.ones_like(data2)+2, data2, c='k', alpha=0.2, linewidths=0)
        plt.plot([2,3], [np.mean(data1), np.mean(data2)], color='k', alpha=0.2)
    MEAN_CONTROL_SALINE_acrossMice=np.mean(control_plot[0])
    MEAN_CONTROL_AGONIST_acrossMice=np.mean(control_plot[1])
    SEM_CONTROL_SALINE_acrossMice=np.std(control_plot[0])/np.sqrt(len(control_plot[0]))
    SEM_CONTROL_AGONIST_acrossMice=np.std(control_plot[1])/np.sqrt(len(control_plot[1]))
    
    plt.plot([2,3], [MEAN_CONTROL_SALINE_acrossMice, MEAN_CONTROL_AGONIST_acrossMice], color='red', linewidth=4)
    plt.vlines(2,MEAN_CONTROL_SALINE_acrossMice-SEM_CONTROL_SALINE_acrossMice, MEAN_CONTROL_SALINE_acrossMice+SEM_CONTROL_SALINE_acrossMice, color='red', linewidths=4)
    plt.vlines(3,MEAN_CONTROL_AGONIST_acrossMice-SEM_CONTROL_AGONIST_acrossMice, MEAN_CONTROL_AGONIST_acrossMice+SEM_CONTROL_AGONIST_acrossMice, color='red', linewidths=4)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylabel('CR rate', size=20)
    plt.xticks([0,1,2,3], size=14)
    plt.xlim(-0.5,3.5)
    # plt.yticks([0,0.025,0.05,0.075,0.1],size=14)
    #plt.yticks([0,0.05,0.1,0.15,0.2],size=14)
    plt.yticks([0,0.2, 0.4, 0.6, 0.8],size=14)
    plt.title('DREADD                                YFP \nSaline              Agonist21             Saline         Agonist21')     
    
    
    
    
    
    
    Anova_df=pd.DataFrame(columns=['feature','virus','injection', 'subject'])
    counter=0
    for a,b,c in zip(dreadd_plot[0],dreadd_plot[1], dreadd_plot[2]):
          Anova_df.at[counter,'feature']=a
          Anova_df.at[counter,'virus']=1
          Anova_df.at[counter,'injection']=0
          Anova_df.at[counter,'subject']=c
          counter+=1
          Anova_df.at[counter,'feature']=b
          Anova_df.at[counter,'virus']=1
          Anova_df.at[counter,'injection']=1
          Anova_df.at[counter,'subject']=c
          counter+=1
    for a,b,c in zip(control_plot[0],control_plot[1], control_plot[2]):
          Anova_df.at[counter,'feature']=a
          Anova_df.at[counter,'virus']=0
          Anova_df.at[counter,'injection']=0
          Anova_df.at[counter,'subject']=c
          counter+=1
          Anova_df.at[counter,'feature']=b
          Anova_df.at[counter,'virus']=0
          Anova_df.at[counter,'injection']=1
          Anova_df.at[counter,'subject']=c
          counter+=1
    
    # #degrees of freedom
    # N = len(Anova_df.feature)
    # df_a = len(Anova_df.virus.unique()) - 1
    # df_b = len(Anova_df.injection.unique()) - 1
    # df_axb = df_a*df_b 
    # df_w = N - (len(Anova_df.virus.unique())*len(Anova_df.injection.unique()))
    
    
    # #Sum of squares
    # grand_mean = Anova_df['feature'].mean()
    # ssq_a = sum([(Anova_df[Anova_df.virus ==l].feature.mean()-grand_mean)**2 for l in Anova_df.virus])
    # ssq_b = sum([(Anova_df[Anova_df.injection ==l].feature.mean()-grand_mean)**2 for l in Anova_df.injection])
    # ssq_t = sum((Anova_df.feature - grand_mean)**2)
    
    # #Sum of squares within
    # dreadd = Anova_df[Anova_df.virus == 1]
    # control = Anova_df[Anova_df.virus == 0]
    # dreadd_injection_means = [dreadd[dreadd.injection == d].feature.mean() for d in dreadd.injection]
    # control_injection_means = [control[control.injection == d].feature.mean() for d in control.injection]
    # ssq_w = sum((control.feature - control_injection_means)**2) +sum((dreadd.feature - dreadd_injection_means)**2)
    
    # #Sum of squares interaction
    # ssq_axb = ssq_t-ssq_a-ssq_b-ssq_w
    
    # #Mean Square A
    # ms_a = ssq_a/df_a
    # #Mean Square B
    # ms_b = ssq_b/df_b
    # #Mean Square AxB
    # ms_axb = ssq_axb/df_axb
    # #Mean Square Within/Error/Residual
    # ms_w = ssq_w/df_w
    
    # #F-ratio
    # f_a = ms_a/ms_w
    # f_b = ms_b/ms_w
    # f_axb = ms_axb/ms_w
    
    # #p-values
    # p_a = stats.f.sf(f_a, df_a, df_w)
    # p_b = stats.f.sf(f_b, df_b, df_w)
    # p_axb = stats.f.sf(f_axb, df_axb, df_w)
    
    
    # #RESULTS
    # results = {'sum_sq':[ssq_a, ssq_b, ssq_axb, ssq_w],
    #            'df':[df_a, df_b, df_axb, df_w],
    #            'F':[f_a, f_b, f_axb, 'NaN'],
    #             'PR(&gt;F)':[p_a, p_b, p_axb, 'NaN']}
    # columns=['sum_sq', 'df', 'F', 'PR(&gt;F)']
    # aov_table1 = pd.DataFrame(results, columns=columns,
    #                           index=['virus', 'injection', 
    #                           'virus:injection', 'Residual'])
    
    # def eta_squared(aov):
    #     aov['eta_sq'] = 'NaN'
    #     aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    #     return aov
    # def omega_squared(aov):
    #     mse = aov['sum_sq'][-1]/aov['df'][-1]
    #     aov['omega_sq'] = 'NaN'
    #     aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    #     return aov
    # eta_squared(aov_table1)
    # omega_squared(aov_table1)
    # print(aov_table1)
    
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    from statsmodels.stats.anova import AnovaRM
    from statsmodels.regression.mixed_linear_model import MixedLMResults
    
    # #add a 'subject_id' to the df
    counter=0
    k=0
    for i in Anova_df.index:
        Anova_df.at[i,'subject']=k
        counter+=1
        if counter%2==0:
            k+=1
    # for each in Anova_df.index:
    #     Anova_df.loc[each,'feature']=np.round(Anova_df.loc[each,'feature'], 5)
    import pingouin as pg
    # Compute the two-way mixed-design ANOVA
    Anova_df['feature'] = pd.to_numeric(Anova_df['feature'])
    aov = pg.mixed_anova(dv='feature', within='injection', between='virus', subject='subject', data=Anova_df)
    # Pretty printing of ANOVA summary
    pg.print_table(aov)
    
    
    
    
    s,p_dreaddSALINE_dreaddAGONIST=stats.ttest_rel(dreadd_plot[0], dreadd_plot[1])
    s,p_controlSALINE_controlAGONIST=stats.ttest_rel(control_plot[0], control_plot[1])
    s,p_dreaddSALINE_controlSALINE=stats.ttest_ind(dreadd_plot[0], control_plot[0])
    s,p_dreaddAGONIST_controlAGONIST=stats.ttest_ind(dreadd_plot[1], control_plot[1])
    
    return p_dreaddSALINE_dreaddAGONIST, p_controlSALINE_controlAGONIST, p_dreaddSALINE_controlSALINE,p_dreaddAGONIST_controlAGONIST