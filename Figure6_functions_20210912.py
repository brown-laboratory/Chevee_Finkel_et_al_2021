###############################################################################
# Figure 6
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
    # - LickPrefCutoff()
    # - example_baseline()
    # - Intention_index()


###############################################################################
# LickPrefCutoff
###############################################################################

def LickPrefCutoff(mice, master_log):
    LickPref=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    LickPref.append(unit_log['LickPref'].values[0])
                    print(unit)
    LickPref=[x  for x in LickPref if not (math.isnan(np.sum(x))) ]

    plt.hist(LickPref, bins=50, color='grey', alpha=0.3, linewidth=0)
    plt.vlines(-0.1, 0,60)
    plt.vlines(0.1, 0,60)
    plt.xlabel('LickPref')
    plt.ylabel('Count')
    
    Right_pref=len([x for x in LickPref if x>0.1])
    Left_pref=len([x for x in LickPref if x<-0.1])
    
    return Right_pref, Left_pref

###############################################################################
# example_baseline
###############################################################################
    
def example_baseline(master_log, example): #'Claustrum31_20191107_TT7clst1'   'Claustrum31_20191114_TT3clst1' 
    unit_log=master_log.loc[np.equal(master_log['unit_name'],example)]

    fig, ax = plt.subplots(1,1)
    counter=0
    for trial in unit_log.index:
        trial_log=unit_log.loc[trial,:]
        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
        if not sum(maskR+maskL):
            spikes=trial_log['StimALigned_spike_times'][(trial_log['StimALigned_spike_times']<0) & (trial_log['StimALigned_spike_times']>-1)]
            plt.scatter(spikes, np.zeros_like(spikes)+counter, c='k', s=0.5) 
            if trial_log['block_type']=='Whisker':
                plt.scatter([0.01,0.02], np.zeros((2))+counter, c='slateblue', s=0.5) 
                counter+=1
            if trial_log['block_type']=='Visual':
                plt.scatter([0.01,0.02], np.zeros((2))+counter, c='tomato', s=0.5) 
                counter+=1
 
        
    plt.xticks([-1, -0.5, 0])
    plt.vlines(0,1,counter-1)
    plt.xlim(-1,0.04)
    plt.ylim(0,counter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylabel('Trials', size=18)
    plt.xlabel('Time (seconds)', size=18)
    plt.title(unit_log['unit_name'].values[0])
    return

###############################################################################
# Intention_index
###############################################################################
    
def Intention_index(master_log, mice):
    from scipy.stats import gaussian_kde

    ORIGINAL=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37']
    REVERSED=['Claustrum18','Claustrum23','Claustrum25'] 
    temp_units1=[] #RightPref, Right trial
    temp_units2=[] #RightPref, Left trials
    temp_units3=[] # LeftPref, Right trials
    temp_units4=[] # LeftPref, Left trials
    names_Right=[]
    names_Left=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log[day_log['LickPref']>0.1]['cluster_name']): #Right
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                print(unit)
                if mouse in ORIGINAL:
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Whisker'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units1.append(np.mean(temp, axis=0))
                    
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Visual'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units2.append(np.mean(temp, axis=0))
                        
                elif mouse in REVERSED:
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Whisker'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units2.append(np.mean(temp, axis=0))
                    
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Visual'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units1.append(np.mean(temp, axis=0))
                names_Right.append(unit_log['unit_name'].values[0])
                        
            for unit in np.unique(day_log[day_log['LickPref']<-0.1]['cluster_name']): #Right
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                if mouse in ORIGINAL:
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Whisker'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units3.append(np.mean(temp, axis=0))
                    
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Visual'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units4.append(np.mean(temp, axis=0))
                        
                elif mouse in REVERSED:
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Whisker'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units4.append(np.mean(temp, axis=0))
                    
                    temp=[]
                    for trial in unit_log[unit_log['block_type']=='Visual'].index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if not sum(maskR+maskL):
                            temp.append(trial_log['-1to3sec_25msecbins_StimAligned'][:41])
                    temp_units3.append(np.mean(temp, axis=0))
                names_Left.append(unit_log['unit_name'].values[0])
                    
    temp_units1=[x  for x,y in zip(temp_units1, temp_units2) if not (math.isnan(np.sum(x)) | math.isnan(np.sum(y)))]
    temp_units2=[y  for x,y in zip(temp_units1, temp_units2) if not (math.isnan(np.sum(x)) | math.isnan(np.sum(y)))]
    temp_units3=[x  for x,y in zip(temp_units3, temp_units4) if not (math.isnan(np.sum(x)) | math.isnan(np.sum(y)))]
    temp_units4=[y  for x,y in zip(temp_units3, temp_units4) if not (math.isnan(np.sum(x)) | math.isnan(np.sum(y)))]
    
    X=[(np.mean(a)-np.mean(b))/(np.mean(a)+np.mean(b)) for a,b in zip(temp_units1, temp_units2)]
    Y=[(np.mean(a)-np.mean(b))/(np.mean(a)+np.mean(b)) for a,b in zip(temp_units3, temp_units4)]
    
    fig,ax=plt.subplots(3,1,figsize=(10,10))

    plt.sca(ax[0])
    density = gaussian_kde(X)
    xs = np.linspace(-1,1,200)
    plt.vlines(0,0,2, linestyles='dotted')
    ax[0].fill_between(xs, 0,density(xs), color='#A1007D', alpha=0.5)
    plt.xlim(-1.05,1.05)
    s,p1=stats.wilcoxon(X,np.zeros_like(X))
    s,p2=stats.wilcoxon(Y,np.zeros_like(Y))
    plt.title('Intention index, RightPref n='+str(len(X))+', p='+str(p1)+'\n LeftPref n='+str(len(Y))+', p='+str(p2))
    
    plt.sca(ax[2])
    y=np.random.random(len(X))
    plt.scatter(X, y, c='#A1007D',s=50, alpha=0.5, linewidths=0)
    plt.vlines(0,0,1, linestyles='dotted')
    plt.xlim(-1.05,1.05)
    
    #Add y
    plt.sca(ax[0])
    density = gaussian_kde(Y)
    xs = np.linspace(-1,1,200)
    plt.vlines(0,0,2, linestyles='dotted')
    ax[0].fill_between(xs, 0,density(xs), color='#00D88A', alpha=0.5)
    plt.xlim(-1.05,1.05)
    plt.sca(ax[1])
    
    box=plt.boxplot([ X,Y], vert=False, widths=0.2,sym='.',medianprops=dict(linestyle='-', linewidth=2.5, color='k'), patch_artist=True)
    for patch, color in zip(box['boxes'], ['#A1007D','#00D88A']):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
    plt.vlines(0,0,4, linestyles='dotted')
    plt.xlim(-1.05,1.05)
    
    plt.sca(ax[2])
    
    y=np.random.random(len(Y))
    plt.scatter(Y, y, c='#00D88A', alpha=0.5, s=50, linewidths=0)
    plt.vlines(0,0,1, linestyles='dotted')
    plt.xlim(-1.05,1.05)
    
    #get the Intention index of example
    example1='Claustrum31_20191114_TT3clst1'
    ExampleIndex=[i for i,x in zip(X,names_Right) if x==example1]
    plt.scatter(ExampleIndex, 0.56, c='r', alpha=1, s=100, linewidths=0)
    
    example2='Claustrum31_20191107_TT7clst1'
    ExampleIndex=[i for i,x in zip(Y,names_Left) if x==example2]
    plt.scatter(ExampleIndex, 0.56, c='r', alpha=1, s=100, linewidths=0)
    
    Right_mean=np.mean(X)
    Right_sem=np.std(X)/np.sqrt(len(X))
    Left_mean=np.mean(Y)
    Left_sem=np.std(Y)/np.sqrt(len(Y))
    s,p1=stats.wilcoxon(X, np.zeros_like(X))
    s,p2=stats.wilcoxon(Y,np.zeros_like(Y))
    s,p3=stats.mannwhitneyu(X,Y)
    
    return Right_mean, Right_sem, Left_mean, Left_sem, p1, p2, p3

###############################################################################
# DREADD_anova_wrapper
###############################################################################

def DREADD_anova_wrapper(mice, master_DREADD, test):   
    
    if test=="Hit_rate":
        def test_function(day_log): #dprime
            from scipy.stats import norm
            HIT_relevant_trials=0
            FA_relevant_trials=0
            HITs=0
            FAs=0
            counter=0
            for trialtype, BlockType,StimType in zip([['SomHit', 'SomFA1','SomFA2'] , ['VisHit', 'VisFA1', 'VisFA2']], ['Whisker','Visual'], ['Stim_Som_NoCue','Stim_Vis_NoCue']):
                block_log=day_log[day_log['block_type']==BlockType]
                HIT_relevant_trials+=len(block_log[block_log['trial_type']==StimType])
                FA_relevant_trials+=len(block_log[block_log['trial_type']!=StimType])
                HITs+=len(block_log[block_log['Stim/Block/Response']==trialtype[0]])
                FAs+=len(block_log[block_log['Stim/Block/Response']==trialtype[1]]) + len(block_log[block_log['Stim/Block/Response']==trialtype[2]])
            Hit_rate=HITs/HIT_relevant_trials
            FA_rate=FAs/FA_relevant_trials
            if (Hit_rate!=0) & (FA_rate!=0):
                d_prime=norm.ppf(Hit_rate)-norm.ppf(FA_rate)
            else:
                d_prime=0
            return Hit_rate
    elif test=='FA_rate_impulsivity':
        def test_function(day_log):
            from scipy.stats import norm
            HIT_relevant_trials=0
            FA_relevant_trials=0
            HITs=0
            FAs=0
            counter=0
            for trialtype, BlockType,StimType in zip(['SomFA2' , 'VisFA2'], ['Whisker','Visual'], ['Stim_Vis_NoCue','Stim_Som_NoCue']):
                block_log=day_log[day_log['block_type']==BlockType]
                FA_relevant_trials+=len(block_log[block_log['trial_type']!=StimType])
                FAs+=len(block_log[block_log['Stim/Block/Response']==trialtype]) 
            FA_rate=FAs/FA_relevant_trials
            return  FA_rate
    elif test=='FA_rate_rule_exclusive':
        def test_function(day_log):
            from scipy.stats import norm
            HIT_relevant_trials=0
            FA_relevant_trials=0
            HITs=0
            FAs=0
            counter=0
            for trialtype, BlockType,StimType in zip(['SomFA3'  , 'VisFA3'], ['Whisker','Visual'], ['Stim_Som_NoCue','Stim_Vis_NoCue']):
                block_log=day_log[day_log['block_type']==BlockType]
                FA_relevant_trials+=len(block_log[block_log['trial_type']!=StimType])
                FAs+=len(block_log[block_log['Stim/Block/Response']==trialtype]) 
            FA_rate=FAs/FA_relevant_trials
            return  FA_rate
    ## Run with Hit rate
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
                       # print(mouse)
                        for day in np.unique(mouse_log['date'])[:8]:
                            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                            #print(day)
                            Rate= test_function(day_log)
                            #print(Hit_rate)
                            if day_log['Injection'].values[0]=='SALINE':
                                  SALINE_data.append([Rate])
                            if day_log['Injection'].values[0]=='AGONIST21':
                                  AGONIST21_data.append([Rate])
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
               # print('_pairedSessions_'+str(p))
                
                XX.append(np.mean(SALINE_data))
                YY.append(np.mean(AGONIST21_data))
                ZZ.append(mouse)
          s,p=stats.wilcoxon(XX, YY)
          #print('pairedMice_'+str(p))
                      
          #Summary_data.append([All_Saline, All_Agonist21, All_mouse])
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


