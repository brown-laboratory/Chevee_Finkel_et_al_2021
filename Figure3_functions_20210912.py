###############################################################################
# Figure 3
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
    # - LDA_figure()
    # - AUC_OptoTag_Modality()
    
##############################################################################
# LDA_figure()
##############################################################################

def LDA_figure(master_log, mice):
    import random
    from sklearn.model_selection import train_test_split
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    from matplotlib.pyplot import figure
    from sklearn.metrics import confusion_matrix
    from scipy.stats import chi2_contingency
    
    
    All_Matrices=[]
    All_TrialTypes=[]
    Mouse_day=[]
    Session_size=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            Day_mod_neurons=0
            print(day)
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            Number_of_units=len(np.unique(day_log['cluster_name']))
            Session_size.append(Number_of_units)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'],np.unique(day_log['cluster_name'])[0][0])]
            Number_of_trials=len(unit_log)
            Matrix=np.zeros((Number_of_trials, Number_of_units))
            unit_count=0
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
    
                TrialTypes=[]
                trial_count=0
                for trial in unit_log.index:
                    trial_log=unit_log.loc[trial,:]
                    data=np.sum(trial_log['-1to3sec_25msecbins_StimAligned'][40:80])
                    Matrix[trial_count, unit_count]=data
                    TrialTypes.append(trial_log['Stim/Block/Response'])
                    trial_count+=1
                unit_count+=1
            All_Matrices.append(Matrix)
            All_TrialTypes.append(TrialTypes)
                
    
    All_Cat_scores=[]
    Category_type=['Lick' ,'Stim', 'Block']
    for Ctype in Category_type:
        All_Scores=[]
        for each,X in enumerate(All_Matrices):#rows: neurons, columns: trials
            if Ctype=='Lick':
                Trial_Category=[]
                for trial in All_TrialTypes[each]:
                    if (trial=='SomHit') | (trial=='VisHit') | (trial=='SomFA') | (trial=='VisFA'):
                        Trial_Category.append(0)
                    elif (trial=='SomMiss') | (trial=='VisMiss') | (trial=='SomCR') | (trial=='VisCR'):
                        Trial_Category.append(1)
            if Ctype=='Stim':
                Trial_Category=[]
                for trial in All_TrialTypes[each]:
                    if (trial=='SomHit') | (trial=='VisCR') | (trial=='SomMiss') | (trial=='VisFA'):
                        Trial_Category.append(0)
                    elif (trial=='SomCR') | (trial=='VisMiss') | (trial=='SomFA') | (trial=='VisHit'):
                        Trial_Category.append(1)
            if Ctype=='Block':
                Trial_Category=[]
                for trial in All_TrialTypes[each]:
                    if (trial=='SomHit') | (trial=='SomCR') | (trial=='SomMiss') | (trial=='SomFA'):
                        Trial_Category.append(0)
                    elif (trial=='VisMiss') | (trial=='VisHit') | (trial=='VisCR') | (trial=='VisFA'):
                        Trial_Category.append(1)
                    
                    
            
            SessionScores=[]
            for i in np.arange(10):
                #random.shuffle(y)
                X_train, X_test, y_train, y_test = train_test_split(
                    X, Trial_Category,
                    test_size=0.2, 
                    stratify=Trial_Category # this makes sure that our training and testing sets both have all classes in y
                )
            
                classifier = LinearDiscriminantAnalysis(n_components=1)
                scores = cross_val_score(classifier, X_train, y_train, cv=5)
                SessionScores.append(scores)
            All_Scores.append(np.mean(SessionScores, axis=0))
            print(each)
        All_Cat_scores.append(All_Scores)
        #sns.distplot(np.mean(All_Scores, axis=1))
    #NOTE: the way it is done here is actually overkill as I do both use only 80% of the data and use cross-validation (80%again)
    
    each=3
    X=All_Matrices[each]
    #plot example session LICK
    Trial_Category=[]
    for trial in All_TrialTypes[each]:
        if (trial=='SomHit') | (trial=='VisHit') | (trial=='SomFA') | (trial=='VisFA'):
            Trial_Category.append(0)
        elif (trial=='SomMiss') | (trial=='VisMiss') | (trial=='SomCR') | (trial=='VisCR'):
            Trial_Category.append(1)
    classifier.fit(X, Trial_Category)
    plt.figure(figsize=(5,5))        
    group1=classifier.transform([x for x,i in zip(X,Trial_Category) if i==0])
    group2=classifier.transform([x for x,i in zip(X,Trial_Category) if i==1])
    plt.hist(group1, color='k', bins=50, range=[-4.5,4.5], alpha=0.8)
    plt.hist(group2, color='lightgrey', bins=50, range=[-4.5,4.5], alpha=0.9)
    plt.xlabel('LDA 1')
    plt.ylabel('Counts')
    plt.title('Lick vs NoLick, n= ' + str(np.size(X,1)) + ' neurons')
    
    #plot example session STIM
    Trial_Category=[]
    for trial in All_TrialTypes[each]:
        if (trial=='SomHit') | (trial=='VisCR') | (trial=='SomMiss') | (trial=='VisFA'):
            Trial_Category.append(0)
        elif (trial=='SomCR') | (trial=='VisMiss') | (trial=='SomFA') | (trial=='VisHit'):
            Trial_Category.append(1)
    classifier.fit(X, Trial_Category)
    plt.figure(figsize=(5,5))        
    group1=classifier.transform([x for x,i in zip(X,Trial_Category) if i==0])
    group2=classifier.transform([x for x,i in zip(X,Trial_Category) if i==1])
    plt.hist(group1, color='cornflowerblue', bins=50, range=[-4.5,4.5], alpha=0.5)
    plt.hist(group2, color='orange', bins=50, range=[-4.5,4.5], alpha=0.5)
    plt.xlabel('LDA 1')
    plt.ylabel('Counts')
    plt.title('TouchStim vs VisualStim, n= ' + str(np.size(X,1)) + ' neurons')
    
    
    #plot example session BLOCK
    Trial_Category=[]
    for trial in All_TrialTypes[each]:
        if (trial=='SomHit') | (trial=='SomCR') | (trial=='SomMiss') | (trial=='SomFA'):
            Trial_Category.append(0)
        elif (trial=='VisMiss') | (trial=='VisHit') | (trial=='VisCR') | (trial=='VisFA'):
            Trial_Category.append(1)
    classifier.fit(X, Trial_Category)
    plt.figure(figsize=(5,5))        
    group1=classifier.transform([x for x,i in zip(X,Trial_Category) if i==0])
    group2=classifier.transform([x for x,i in zip(X,Trial_Category) if i==1])
    plt.hist(group1, color='slateblue', bins=50, range=[-4.5,4.5], alpha=0.5)
    plt.hist(group2, color='tomato', bins=50, range=[-4.5,4.5], alpha=0.5)
    plt.xlabel('LDA 1')
    plt.ylabel('Counts')
    plt.title('TouchBlock vs VisualBlock, n= ' + str(np.size(X,1)) + ' neurons')
    
            
    #Summary plot
    fig, ax=plt.subplots(1,1, figsize=(4,8))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.sca(ax)
    for a,b,c in zip(All_Cat_scores[1], All_Cat_scores[2], All_Cat_scores[0]):
        plt.plot([1,2,3], [np.mean(a),np.mean(b),np.mean(c)], color='k', alpha=0.1)
    plt.plot([1,2,3], np.mean(np.mean([All_Cat_scores[1], All_Cat_scores[2], All_Cat_scores[0]],axis= 1), axis=1), color='k', linewidth=2)
    plt.vlines(1,np.mean(All_Cat_scores[1])-np.std(np.mean(All_Cat_scores[1], axis=1))/np.sqrt(len(All_Cat_scores[1])), np.mean(All_Cat_scores[1])+np.std(np.mean(All_Cat_scores[1],axis=1))/np.sqrt(len(All_Cat_scores[1])), linewidth=2)
    plt.vlines(2,np.mean(All_Cat_scores[2])-np.std(np.mean(All_Cat_scores[2], axis=1))/np.sqrt(len(All_Cat_scores[2])), np.mean(All_Cat_scores[2])+np.std(np.mean(All_Cat_scores[2],axis=1))/np.sqrt(len(All_Cat_scores[2])), linewidth=2)
    plt.vlines(3,np.mean(All_Cat_scores[0])-np.std(np.mean(All_Cat_scores[0], axis=1))/np.sqrt(len(All_Cat_scores[0])), np.mean(All_Cat_scores[0])+np.std(np.mean(All_Cat_scores[0],axis=1))/np.sqrt(len(All_Cat_scores[0])), linewidth=2)
    plt.ylim(0.5,1)
    plt.ylabel('Mean classification accuracy')
    s,p_anova=sp.stats.f_oneway(np.mean(All_Cat_scores[0], axis=1), np.mean(All_Cat_scores[1], axis=1), np.mean(All_Cat_scores[2], axis=1))
    s,p_lick_vs_stim=sp.stats.wilcoxon(np.mean(All_Cat_scores[0], axis=1) , np.mean(All_Cat_scores[1], axis=1))
    s,p_lick_vs_block=sp.stats.wilcoxon(np.mean(All_Cat_scores[0], axis=1) , np.mean(All_Cat_scores[2], axis=1))
    s,p_stim_vs_block=sp.stats.wilcoxon(np.mean(All_Cat_scores[1], axis=1) , np.mean(All_Cat_scores[2], axis=1))
    
    #Summary data
    meanSessionSize=np.mean(Session_size)
    minSessionSize= np.min( Session_size)
    maxSessionSize=np.max( Session_size)
    semSessionSize=np.std( Session_size)/np.sqrt(len( Session_size))
    
    #lick
    meanScoreLick=np.mean(np.mean(All_Cat_scores[0], axis=1))
    semScoreLick=np.std(np.mean(All_Cat_scores[0], axis=1))/np.sqrt(len(np.mean(All_Cat_scores[0], axis=1)))
    
    #stim
    meanScoreStim= np.mean(np.mean(All_Cat_scores[1], axis=1))
    semScoreStim= np.std(np.mean(All_Cat_scores[1], axis=1))/np.sqrt(len(np.mean(All_Cat_scores[1], axis=1)))
    
    #block
    meanScoreBlock=np.mean(np.mean(All_Cat_scores[2], axis=1))
    semScoreBlock=np.std(np.mean(All_Cat_scores[2], axis=1))/np.sqrt(len(np.mean(All_Cat_scores[2], axis=1)))
    
    return  meanSessionSize, minSessionSize,  maxSessionSize,  semSessionSize, meanScoreLick, semScoreLick, meanScoreStim, semScoreStim, meanScoreBlock, semScoreBlock, p_anova, p_lick_vs_stim, p_lick_vs_block, p_stim_vs_block


##############################################################################
# AUC_mean and AUC_scatter OptoTag Only
##############################################################################
def AUC_OptoTag_Modality(mice, master_log, BinNumber):
    
    
    temp_units1=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']=='OptoTag']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_units1.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
    
    temp_units2=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']=='OptoTag']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    temp_units2.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
    
    plt.style.use('default')
    fig, ax=plt.subplots(1,1)                
    SEM1=sp.stats.sem(temp_units1, axis=0)
    SomHit, =ax.plot(np.arange(159), np.mean(temp_units1, axis=0), color='cornflowerblue')
    ax.fill_between(np.arange(159),np.mean(temp_units1, axis=0)+SEM1, np.mean(temp_units1, axis=0)-SEM1, color='cornflowerblue', alpha=0.1)
                      
    SEM2=sp.stats.sem(temp_units2, axis=0)
    VisHit, =ax.plot(np.arange(159), np.mean(temp_units2, axis=0), color='orange')
    ax.fill_between(np.arange(159),np.mean(temp_units2, axis=0)+SEM1, np.mean(temp_units2, axis=0)-SEM1, color='orange', alpha=0.1)
      
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
    ax.legend([SomHit, VisHit],['Touch','Visual'])
    fig.suptitle('Mean AUC - Stim Aligned - Hit vs CR\n All units, n=' + str(len(temp_units1)))
    
    
    #Scatter meanAUC (500msec after onset) -  optotag
    
    temp_units3=[] 
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']=='OptoTag']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]
                        temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units3.append(0.5)
    temp_units4=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']=='OptoTag']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0]
                        temp_units4.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units4.append(0.5)
                         
    plt.style.use('default')
    fig, ax=plt.subplots(1,1, figsize=(8,8))
    ax.fill([0.3,1,1],[0.3,1,0.3], 'cornflowerblue', alpha=0.3)
    ax.fill([0.3,0.3,1],[0.3,1,1], 'orange', alpha=0.2)
    ax.set_xlim(0.3,1)
    ax.set_ylim(0.3,1)
    ax.scatter(temp_units3, temp_units4, c='k',s=30, linewidths=0)
    ax.plot([0.3,1],[0.3,1], color='k')
    ax.vlines(0.5, 0.3,1, linestyles='dotted')
    ax.hlines(0.5,0.3,1, linestyles='dotted')
    ax.set_xlabel('Mean AUC Touch')
    ax.set_ylabel('Mean AUC Visual')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    s,p=sp.stats.wilcoxon(temp_units3, temp_units4)
    s,p1=sp.stats.wilcoxon([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))] , [y for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])

    fig.suptitle('Mean AUC ('+str(BinNumber*25)+'msec) - Stim Aligned - Hit vs CR \n OptoTag, n='+str(len(temp_units3))+ ' wilcoxon p='+str(p)[0:5] +'\n Sig: n='+
                 str(len([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])) +' wilcoxon p=' + str(p1)[0:5])
    
    
    #numers
    Touch_mean=np.mean(temp_units3)
    Touch_std=np.std(temp_units3)
    Touch_sem=np.std(temp_units3)/np.sqrt(len(temp_units3))
    Visual_mean=np.mean(temp_units4)
    Visual_std=np.std(temp_units4)
    Visual_sem=np.std(temp_units4)/np.sqrt(len(temp_units4))
    
    plt.plot(Touch_mean,  Visual_mean, color='k', marker='o')
    plt.vlines(Touch_mean,  Visual_mean-Visual_sem,  Visual_mean+Visual_sem)
    plt.hlines(Visual_mean, Touch_mean-Touch_sem, Touch_mean+Touch_sem)
    return  Touch_mean, Touch_sem, Visual_mean, Visual_sem
