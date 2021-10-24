###############################################################################
# Figure Reviewers
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
    # - 
    
    
###############################################################################
# Figure S7A,B
###############################################################################
def Tetrode_map_LickPref(mice, master_log):    
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        plt.figure()
        grid = plt.GridSpec(len(np.unique(mouse_log['date'])),8, wspace=0, hspace=0)
        for i,day in enumerate(np.unique(mouse_log['date'])):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            TT1=[]
            TT2=[]
            TT3=[]
            TT4=[]
            TT5=[]
            TT6=[]
            TT7=[]
            TT8=[]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                
                if unit_log['cluster_name'].values[0][0][:3]=='TT1':
                    TT1.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT2':
                    TT2.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT3':
                    TT3.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT4':
                    TT4.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT5':
                    TT5.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT6':
                    TT6.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT7':
                    TT7.append(unit_log['LickPref'].values[0])
                elif unit_log['cluster_name'].values[0][0][:3]=='TT8':
                    TT8.append(unit_log['LickPref'].values[0])
            for j,TT in enumerate ([TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8]):
                ax=plt.subplot(grid[i,j])
                RightPref=len([x for x in TT if x>0.1])
                LeftPref=len([x for x in TT if x<-0.1])
                NoPref=len([x for x in TT if ((x<0.1) &(x>-0.1))])
                ax.pie([NoPref, RightPref,LeftPref], 
                    colors=['darkgrey','#A1007D', '#00D88A'], 
                    startangle=90 , 
                    counterclock=False )
        plt.suptitle(mouse)
    #    plt.savefig('Tetrode map LickPref ' + mouse + '.pdf')
    #    plt.close()
    return


###############################################################################
# Does FIgure 4C depend on 10 outlyers?
###############################################################################
    
def Intention_noOutlyers(mice, master_log):
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
    

    # Both distributions are normal, so we will use Zscore as a means to exclude outliers
    sp.stats.normaltest(X)
    sp.stats.normaltest(Y)
    
    names_Right=[names_Right[i] for i in np.argsort(X)]
    names_Left=[names_Left[i] for i in np.argsort(Y)]
    X=np.sort(X)
    Y=np.sort(Y)
    
    Zscored_X=stats.zscore(X)
    Zscored_Y=stats.zscore(Y)
    fig, ax =plt.subplots(1,1,figsize=(5,5))
    plt.hist(Zscored_X, alpha=0.5)
    plt.hist(Zscored_Y, alpha=0.5)
    plt.legend(['contra', 'ipsi'])
    plt.vlines(-3,0,50)
    plt.vlines(3,0,50)
    #Number of outliers in X:4 ; Y:3
    sum((Zscored_X>3) | (Zscored_X<-3))
    sum((Zscored_Y>3) | (Zscored_Y<-3))
    #remove those outliers
    X=[x for x,y in zip(X, Zscored_X) if ((y>-3) & (y<3))]
    Y=[x for x,y in zip(Y, Zscored_Y) if ((y>-3) & (y<3))]
    names_Right=[x for x,y in zip(names_Right, Zscored_X) if ((y>-3) & (y<3))]
    names_Left=[x for x,y in zip(names_Left, Zscored_Y) if ((y>-3) & (y<3))]
    plt.xlabel('Zscored LII')
    plt.ylabel('Counts')
    plt.title('3 and 4 outliers were removed')
    #Replot: the test should be a ttest, not a wilcoxon
    from scipy.stats import gaussian_kde
    
    fig,ax=plt.subplots(3,1,figsize=(10,10))
    
    plt.sca(ax[0])
    density = gaussian_kde(X)
    xs = np.linspace(-1,1,200)
    plt.vlines(0,0,2, linestyles='dotted')
    ax[0].fill_between(xs, 0,density(xs), color='#A1007D', alpha=0.5)
    plt.xlim(-1.05,1.05)
    s,p1=stats.ttest_rel(X,np.zeros_like(X))
    s,p2=stats.ttest_rel(Y,np.zeros_like(Y))
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
    s,p1=stats.ttest_rel(X, np.zeros_like(X))
    s,p2=stats.ttest_rel(Y,np.zeros_like(Y))
    s,p3=stats.ttest_ind(X,Y)
    #results: Right_mean:0.048, Right_sem:0.00970, Left_mean:-0.0288, Left_sem:0.0108
    
    return Right_mean, Right_sem, Left_mean, Left_sem, p1, p2, p3