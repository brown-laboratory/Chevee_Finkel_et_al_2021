###############################################################################
# Figure S3
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
    # - Hit_vs_FA()
    # - TouchHit_vs_VisualHit_OPTO()
    # - TouchHit_vs_VisualHit()
    # - meanAUC_category()
    # - All_AUC_modality()
    # - Onset_cumfreq()
    # - Onset_scatter_bar()
    
##############################################################################
# Hit_vs_FA()
##############################################################################

def Hit_vs_FA(master_log, mice):    
    temp_units1=[]
    temp_units2=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            print(mouse + day[0])
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                
                temp_units1_Hits=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_Hits.append(len(spikes))
                for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_Hits.append(len(spikes))       
                temp_units1.append(np.mean(temp_units1_Hits))
                
                temp_units1_FAs=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='SomFA'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_FAs.append(len(spikes))
                for trial in unit_log[unit_log['Stim/Block/Response']=='VisFA'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_FAs.append(len(spikes))       
                temp_units2.append(np.mean(temp_units1_FAs))
    
    fig,ax=plt.subplots(1,1,figsize=(8,8))   
    ax.scatter([np.log(x) for x in temp_units1],[np.log(x) for x in temp_units2], s=50, color='k', alpha=0.5, linewidths=0)
    ax.plot([np.log(0.1),np.log(85)],[np.log(0.1),np.log(85)], color='grey')
    ax.fill([np.log(0.1),np.log(85),np.log(85)],[np.log(0.1),np.log(0.1),np.log(85)], 'cornflowerblue', alpha=0.3)
    ax.fill([np.log(0.1),np.log(85),np.log(0.1)],[np.log(0.1),np.log(85),np.log(85)], 'orange', alpha=0.3)
    s,p3=sp.stats.wilcoxon(temp_units1,temp_units2)
    ax.set_xlabel('Mean response during touch associated lick trials (spikes/sec)', size=16)
    ax.set_ylabel('Mean response during vision associated lick trials (spikes/sec)', size=16)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Spontaneous Licks - wilcoxon - p<'+str(p3)[:5])
    plt.sca(ax)
    plt.xticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.yticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.xlim(np.log(0.1),np.log(85))
    plt.ylim(np.log(0.1),np.log(85))
    Mean1=np.mean(temp_units1)
    SEM1=np.std(temp_units1)/np.sqrt(len(temp_units1))
    Mean2=np.mean(temp_units2)
    SEM2=np.std(temp_units2)/np.sqrt(len(temp_units2))
    
    plt.plot(np.log(Mean1),np.log(Mean2), marker='.', color='r')
    plt.vlines(np.log(Mean1), np.log(Mean2-SEM2), np.log(Mean2+SEM2), colors='r')
    plt.plot([np.log(Mean1-SEM1),np.log(Mean1+SEM1)],[np.log(Mean2), np.log(Mean2)], color='r')
    return Mean1, SEM1, Mean2, SEM2, p3
    
##############################################################################
# TouchHit_vs_VisualHit_OPTO()
##############################################################################

def  TouchHit_vs_VisualHit_OPTO(master_log, mice):
    temp_units1=[]
    temp_units2=[]

    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            print(mouse + day[0])
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                if unit_log['Category'].values[0]!='OptoTag':
                    continue
                temp_units1_SomHit=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_SomHit.append(len(spikes))
                temp_units1.append(np.mean(temp_units1_SomHit))
                
                temp_units1_VisHit=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_VisHit.append(len(spikes))
                temp_units2.append(np.mean(temp_units1_VisHit))
                
    
    
    fig,ax=plt.subplots(1,1,figsize=(8,8))   
    ax.scatter([np.log(x) for x in temp_units1],[np.log(x) for x in temp_units2], s=50, color='g', alpha=0.5, linewidths=0)
    ax.plot([np.log(0.1),np.log(85)],[np.log(0.1),np.log(85)], color='grey')
    ax.fill([np.log(0.1),np.log(85),np.log(85)],[np.log(0.1),np.log(0.1),np.log(85)], 'cornflowerblue', alpha=0.3)
    ax.fill([np.log(0.1),np.log(85),np.log(0.1)],[np.log(0.1),np.log(85),np.log(85)], 'orange', alpha=0.3)
    s,p3=sp.stats.wilcoxon(temp_units1,temp_units2)
    ax.set_xlabel('Mean response during touch associated lick trials (spikes/sec)', size=16)
    ax.set_ylabel('Mean response during vision associated lick trials (spikes/sec)', size=16)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Spontaneous Licks - wilcoxon - p<'+str(p3)[:5])
    plt.sca(ax)
    plt.xticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.yticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.xlim(np.log(0.1),np.log(85))
    plt.ylim(np.log(0.1),np.log(85))
    Mean1=np.mean(temp_units1)
    SEM1=np.std(temp_units1)/np.sqrt(len(temp_units1))
    Mean2=np.mean(temp_units2)
    SEM2=np.std(temp_units2)/np.sqrt(len(temp_units2))
    plt.plot(np.log(Mean1),np.log(Mean2), marker='.', color='r')
    plt.vlines(np.log(Mean1), np.log(Mean2-SEM2), np.log(Mean2+SEM2), colors='r')
    plt.plot([np.log(Mean1-SEM1),np.log(Mean1+SEM1)],[np.log(Mean2), np.log(Mean2)], color='r')
    #Results: mean_touch=19.64 , sem_touch=1.728 , mean_vis=20.546 , sem_vis=1.847, n=73
    return Mean1, SEM1, Mean2, SEM2, p3


##############################################################################
# TouchHit_vs_VisualHit()
##############################################################################
    
def  TouchHit_vs_VisualHit(master_log, mice):
    temp_units1=[]
    temp_units2=[]

    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            print(mouse + day[0])
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                
                temp_units1_SomHit=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_SomHit.append(len(spikes))
                temp_units1.append(np.mean(temp_units1_SomHit))
                
                temp_units1_VisHit=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    spikes=[x[0][0][0] for x in trial_log['StimALigned_spike_times'] if ((x>0) & (x<1))]
                    temp_units1_VisHit.append(len(spikes))
                temp_units2.append(np.mean(temp_units1_VisHit))
                
             
    
    fig,ax=plt.subplots(1,1,figsize=(8,8))   
    ax.scatter([np.log(x) for x in temp_units1],[np.log(x) for x in temp_units2], s=50, color='k', alpha=0.5, linewidths=0)
    ax.plot([np.log(0.1),np.log(85)],[np.log(0.1),np.log(85)], color='grey')
    ax.fill([np.log(0.1),np.log(85),np.log(85)],[np.log(0.1),np.log(0.1),np.log(85)], 'cornflowerblue', alpha=0.3)
    ax.fill([np.log(0.1),np.log(85),np.log(0.1)],[np.log(0.1),np.log(85),np.log(85)], 'orange', alpha=0.3)
    s,p3=sp.stats.wilcoxon(temp_units1,temp_units2)
    ax.set_xlabel('Mean response during touch associated lick trials (spikes/sec)', size=16)
    ax.set_ylabel('Mean response during vision associated lick trials (spikes/sec)', size=16)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('Spontaneous Licks - wilcoxon - p<'+str(p3)[:5])
    plt.sca(ax)
    plt.xticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.yticks([np.log(x) for x in [0.1,5,10,20,40,80]], ['0','5','10','20','40','80'], size=14)
    plt.xlim(np.log(0.1),np.log(85))
    plt.ylim(np.log(0.1),np.log(85))
    Mean1=np.mean(temp_units1)
    SEM1=np.std(temp_units1)/np.sqrt(len(temp_units1))
    Mean2=np.mean(temp_units2)
    SEM2=np.std(temp_units2)/np.sqrt(len(temp_units2))
    plt.plot(np.log(Mean1),np.log(Mean2), marker='.', color='r')
    plt.vlines(np.log(Mean1), np.log(Mean2-SEM2), np.log(Mean2+SEM2), colors='r')
    plt.plot([np.log(Mean1-SEM1),np.log(Mean1+SEM1)],[np.log(Mean2), np.log(Mean2)], color='r')
    #Results: mean_touch=9.403 , sem_touch=0.492 , mean_vis=9.356 , sem_vis=0.481, n=545
    return Mean1, SEM1, Mean2, SEM2, p3


##############################################################################
# meanAUC_category()
##############################################################################

def meanAUC_category(master_log, mice, Category, BinNumber):


    temp_units3=[] 
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']==Category]):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]
                        temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units3.append(0.5)
    temp_units4=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name'][day_log['Category']==Category]):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
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
    fig.suptitle('Mean AUC ('+str(BinNumber*25)+'msec) - Stim Aligned - Hit vs CR \n OptoTag, n='+str(len(temp_units3))+ ' wilcoxon p<'+str(p)[0:5])
    
    #numers
    Som_mean=np.mean(temp_units3)
    Som_sem=np.std(temp_units3)/np.sqrt(len(temp_units3))
    Vis_mean=np.mean(temp_units4)
    Vis_sem=np.std(temp_units4)/np.sqrt(len(temp_units4))
    plt.plot(Som_mean, Vis_mean, color='b', marker='o')
    plt.vlines(Som_mean, Vis_mean- Vis_sem, Vis_mean+ Vis_sem) #np.log(Mean_CR-STD_CR)
    plt.hlines(Vis_mean, Som_mean-Som_sem, Som_mean+Som_sem) #np.log(Mean_Hit-STD_Hit)
    
    
    return Som_mean, Som_sem, Vis_mean, Vis_sem, p


##############################################################################
# All_AUC_modality()
##############################################################################
    
def All_AUC_modality(master_log, mice, BinNumber):
    temp_units3=[] 
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    if unit_log['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]==1:
                        Onset_bin=unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0]
                        temp_units3.append(np.mean(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0][Onset_bin:Onset_bin+BinNumber]))
                    else:
                         temp_units3.append(0.5)
    temp_units4=[]
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
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
    fig.suptitle('Mean AUC ('+str(BinNumber*25)+'msec) - Stim Aligned - Hit vs CR \n OptoTag, n='+str(len(temp_units3))+ ' wilcoxon p<'+str(p)[0:5])
    
    #numers
    Som_mean=np.mean(temp_units3)
    Som_sem=np.std(temp_units3)/np.sqrt(len(temp_units3))
    Vis_mean=np.mean(temp_units4)
    Vis_sem=np.std(temp_units4)/np.sqrt(len(temp_units4))
    plt.plot(Som_mean, Vis_mean, color='b', marker='o')
    plt.vlines(Som_mean, Vis_mean- Vis_sem, Vis_mean+ Vis_sem) #np.log(Mean_CR-STD_CR)
    plt.hlines(Vis_mean, Som_mean-Som_sem, Som_mean+Som_sem) #np.log(Mean_Hit-STD_Hit)
    
    
    #numers SIG
    Touch_mean=np.mean([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])
    Touch_sem=np.std([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])/np.sqrt(len([x for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))]))
    Visual_mean=np.mean([y for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])
    Visual_sem=np.std([y for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])/np.sqrt(len([y for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))]))
    s,p=sp.stats.wilcoxon([x for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))],[y for x,y in zip(temp_units3, temp_units4) if ((x!=0.5) & (y!=0.5))])
    Numbers_SIG=[Touch_mean,Touch_sem,Visual_mean,Visual_sem, p]
    #numers Excited
    Touch_mean=np.mean([x for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))])
    Touch_sem=np.std([x for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))])/np.sqrt(len([x for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))]))
    Visual_mean=np.mean([y for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))])
    Visual_sem=np.std([y for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))])/np.sqrt(len([y for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))]))
    s,p=sp.stats.wilcoxon([x for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))] ,[y for x,y in zip(temp_units3, temp_units4) if ((x>0.5) & (y>0.5))])
    Numbers_EXC=[Touch_mean,Touch_sem,Visual_mean,Visual_sem, p]
    
    #numers Inhibited
    Touch_mean=np.mean([x for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))])
    Touch_sem=np.std([x for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))])/np.sqrt(len([x for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))]))
    Visual_mean=np.mean([y for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))])
    Visual_sem=np.std([y for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))])/np.sqrt(len([y for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))]))
    s,p=sp.stats.wilcoxon([x for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))],[y for x,y in zip(temp_units3, temp_units4) if ((x<0.5) & (y<0.5))])
    Numbers_INH=[Touch_mean,Touch_sem,Visual_mean,Visual_sem, p]
    
    return Som_mean, Som_sem, Vis_mean, Vis_sem, p,Numbers_SIG, Numbers_EXC, Numbers_INH

##############################################################################
# Onset_cumfreq()
##############################################################################
    
def Onset_cumfreq(master_log, mice):
    temp_units1=[]
    temp_units2=[]
    temp_median1=[];
    temp_median2=[]
    count=0
    
    for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                        print(unit)
                        temp_units1.append(unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
                        temp_units2.append(unit_log['unit_Onset_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
                        temp_trial1=[]
                           
                        for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                            trial_log=unit_log.loc[trial,:]
                            temp_trial1.append(trial_log['FirstLick']) 
                        temp_median1.append(np.median(temp_trial1))
                        temp_trial2=[]                    
                        for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                            trial_log=unit_log.loc[trial,:]
                            temp_trial2.append(trial_log['FirstLick']) 
                        temp_median2.append(np.median(temp_trial2))
                        count+=1
                        
    
    x=[x for i,x in enumerate(temp_units1) if ((x>40) & (temp_units2[i]>40))] #TouchHit DP onsets
    y=[y for i,y in enumerate(temp_units2) if ((y>40) & (temp_units1[i]>40))] #VisualHit DP onsets
   
    
    fig, ax=plt.subplots(1,1,figsize=(8,8))
    plt.sca(ax)
    Cum=stats.cumfreq(x, numbins=100)
    A=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    SomDP=plt.plot(A,Cum.cumcount/np.size(x), color='cornflowerblue')
    
    Cum=stats.cumfreq(y, numbins=100)
    B=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    VisDP=plt.plot(B,Cum.cumcount/np.size(y), color='orange')
    
    plt.vlines(40, 0,1, color='k', linestyles='dotted')
    plt.xticks([40,80,120,160],['0','1','2', '3'], size=16)
    plt.yticks([0,0.5,1],['0','0.5','1'], size=16)
    plt.legend(['Touch','Visual'])
    plt.xlabel('DP onset (seconds)', size=20)
    plt.ylabel('ECDF', size=20)
    s,p=stats.wilcoxon(x, y)
    plt.title('ORIGINAL, Wilcoxon on neurons with both, p<'+str(p))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    Som_mean=(np.mean(x)-40)*25
    Som_sem=(np.std([i-40 for i in x]))*25/np.sqrt(len(x))
    Vis_mean=(np.mean(y)-40)*25
    Vis_sem=(np.std([i-40 for i in y]))*25/np.sqrt(len(y))
    return Som_mean, Som_sem, Vis_mean, Vis_sem, x, y, temp_units1, temp_units2, temp_median1, temp_median2, p


##############################################################################
# Onset_scatter_bar()
##############################################################################
    
def Onset_scatter_bar(master_log, mice,  x, y, temp_units1, temp_units2, temp_median1, temp_median2):
    x_median=[x for i,x in enumerate(temp_median1) if ((temp_units1[i]>40) & (temp_units2[i]>40))] #TouchHit medianRT for units with non-zero DP onsets for both
    y_median=[x for i,x in enumerate(temp_median2) if ((temp_units1[i]>40) & (temp_units2[i]>40))] #VisualHit medianRT for units with non-zero DP onsets for both
    median_c=[1 if x>y_median[i] else 0 for i,x in enumerate(x_median) ] #1: TouchHitmedianRT>VisualHitmedianRT
    
    x_touchRT=[m for i,m in enumerate(x) if median_c[i]] #subset of TouchHit DP onsets: only TouchHitmedianRT>VisualHitmedianRT
    y_touchRT=[m for i,m in enumerate(y) if median_c[i]] #subset of VisualHit DP onsets: only TouchHitmedianRT>VisualHitmedianRT
    x_visualRT=[m for i,m in enumerate(x) if not median_c[i]] #subset of TouchHit DP onsets: only TouchHitmedianRT<VisualHitmedianRT
    y_visualRT=[m for i,m in enumerate(y) if not median_c[i]] #subset of VisualHit DP onsets: only TouchHitmedianRT<VisualHitmedianRT
        
    #S1-Ba (scatter)
    fig, ax=plt.subplots(1,1,figsize=(8,8))
    Touch=ax.scatter(x_touchRT,y_touchRT, c='cornflowerblue', alpha=0.5, linewidths=0)
    Visual=ax.scatter(x_visualRT,y_visualRT, c='orange', alpha=0.5, linewidths=0)
    ax.legend([Touch, Visual], ['long TouchRT','long VisualRT'])
    ax.plot(np.arange(40,160),np.arange(40,160), c='k')
    ax.set_title('ORIGINAL, n=' + str(len(x_visualRT)+len(x_touchRT)))
    ax.set_xticks([39.5,79.5, 119.5, 159.5])
    ax.set_xticklabels(['0', '1', '2', '3'], size=16)
    ax.set_yticks([39.5,79.5, 119.5, 159.5])
    ax.set_yticklabels(['0', '1', '2', '3'], size=16)
    ax.set_xlabel('DP onset (s) - Touch trials', size=20)
    ax.set_ylabel('DP onset (s) - Visual trials', size=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(39.5,160)
    ax.set_ylim(39.5,160)
    
    #S1-Bb (bar plot)
    fig, ax=plt.subplots(1,1,figsize=(8,8))
    #cos(pi/4)=0.707
    touchRT_h=[(x_touchRT[i]-y_touchRT[i])*0.707 for i,m in enumerate(x_touchRT)]
    visualRT_h=[(x_visualRT[i]-y_visualRT[i])*0.707 for i,m in enumerate(x_visualRT)]
    ax.vlines(0, 0,60, color='k')
    Visual=ax.hist(visualRT_h, color='orange', bins=20, range=(-100,100), alpha=0.5)
    ax.vlines(np.median(visualRT_h), 0,50, color='orange')
    Touch=ax.hist(touchRT_h, color='cornflowerblue', bins=20, range=(-100,100), alpha=0.5)
    ax.vlines(np.median(touchRT_h), 0,50, color='cornflowerblue')
    
    v,p=sp.stats.mannwhitneyu(touchRT_h, visualRT_h)
    ax.set_title('DP onset by RT - ORIGINAL, Utest:'+str(p)[:4])
    ax.set_xlabel('Longer Visual DPonset <-- 0 --> Longer Touch DPonset')
    ax.set_ylabel('Counts')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return p


