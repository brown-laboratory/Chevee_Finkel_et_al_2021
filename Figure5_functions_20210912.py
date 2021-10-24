###############################################################################
# Figure 5
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
    # - Example_XCorr_pair()
    # - XCorr_example()
    # - Claustrum_v_S1_XCorr()
    # - Quandrant_plot()
    # - PiePlot_LickPref()
    
    

    
    
###############################################################################
# Example_XCorr_pair
###############################################################################

def Example_XCorr_pair(master_log, XCorr_df, pair):
    
    Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
    Filter=np.array([0.05,0.25,0.4,0.25,0.05])
    Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
    plt.figure()
    plt.plot(np.arange(1004),Subtracted_XCorr, color='k',alpha=0.5)
    plt.xticks([0,251,502, 753, 1004], ['-1', '-0.5','0','0.5','1'])
    plt.yticks( [0.02,0.04,0.06,0.08], [x*2.5 for x in [0.02,0.04,0.06,0.08]])
    plt.title(XCorr_df.loc[pair,'Neuron1_name'] +'\n'+ XCorr_df.loc[pair,'Neuron2_name'])
    
    
    Name1=XCorr_df.loc[pair, 'Neuron1_name']
    Name2=XCorr_df.loc[pair, 'Neuron2_name']
    
    
            
    day_log=master_log[(master_log['unit_name']== Name1) | (master_log['unit_name']== Name2)]
    
    Trial_numbers=np.unique(day_log['trial_num'])
    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
    Units_number=2
    #create data frame that will be used to store and compile spikes from all pairs and all trials
    day_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
    for row in day_df.index:
        for col in day_df.index:
            day_df.at[row,col]=np.zeros((1,1000))
    Number_of_counted_trials=   len(Trial_numbers)-1 #we will discount each trial we remove to check that we have enough and then normalize by the right amount     
    
    Spikes_NeuronPair=[] #Neuron 1 an Neuron2 are alternating rows
    for m,trial in enumerate(Trial_numbers[:-1]): #here we need to enumerate to grab the next trial when doing the control (instead of shuffle) (go until the second to last because we are taking '+1')
        Units_log=day_log[day_log['trial_num']==trial] #represent one trial, on ow per neuron
    
        for i in Units_log.index:
            spikes1=Units_log.loc[i,'spike_times']
            spikes1=spikes1[spikes1<Units_log.loc[i,'stim_onset']] #limited to pre-stim
            Binary_spikes1=np.zeros(math.ceil(Units_log.loc[i,'stim_onset']*1000))#
            for spike in spikes1:
                Binary_spikes1[int(spike*1000)]=True
            Binary_spikes1= Binary_spikes1[-1000:]# limit to only the last 1sec before stim
            Spikes_NeuronPair.append(Binary_spikes1)
    
    #DATA
    prev_timestamps=[0,0]
    count=0
    fig,ax= plt.subplots(1,1, figsize=(4,8))            
    for j,row in enumerate(Spikes_NeuronPair):
        if j%2==0:
            color='c'
            marker='.'
            count+=4
        else:
            color='m'
            count+=2
            marker='.'
        timestamps=[i for i,x in enumerate(row) if x]
        plt.scatter(timestamps,np.zeros_like(timestamps)+count, c=color, marker=marker , s=2)
        
        if j%2==1:
            for spike in timestamps:
                if spike in prev_timestamps:
                    plt.vlines(spike, count-20, count+20, color='k', linewidth=4)
        
        prev_timestamps=timestamps
    plt.xlim(0,1000)
    plt.ylim(0,count)
    plt.xticks([0,250,500,750,999], ['-1','-0.75','-0.5','-0.25','0'])
    plt.yticks(np.arange(0,count,300), ['0','50','100','150','200','250', '300','350','400','450'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time (secons)')
    plt.ylabel('Trials')
    
    #Zoom
    fig,ax= plt.subplots(1,1, figsize=(4,8))  
    count=0
    prev_timestamps=[0,0]
    for j,row in enumerate(Spikes_NeuronPair):
        if j%2==0:
            color='c'
            marker='.'
            count+=4
        else:
            color='m'
            count+=2
            marker='.'
        timestamps=[i for i,x in enumerate(row) if x]
        plt.scatter(timestamps,np.zeros_like(timestamps)+count, c=color, marker=marker , s=10)
        
        if j%2==1:
            for spike in timestamps:
                if spike in prev_timestamps:
                    plt.vlines(spike, count-2, count, color='k', linewidth=4)
        
        prev_timestamps=timestamps
    plt.xlim(0,1000)
    plt.ylim(210*6,(240*6)+1)
    plt.xticks([0,250,500,750,999], ['-1','-0.75','-0.5','-0.25','0'])
    plt.yticks(np.arange(1260,1441,60), ['210','220','230','240'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time (secons)')
    plt.ylabel('Trials')
    
    
    
    #SHIFT control
    prev_timestamps=[0,0]
    count=0
    fig,ax= plt.subplots(1,1, figsize=(4,8))            
    for j,row in enumerate(Spikes_NeuronPair):
        if j%2==0:
            color='c'
            marker='.'
            count+=4
        else:
            color='m'
            count+=2
            marker='.'
        timestamps=[i for i,x in enumerate(row) if x]
        plt.scatter(timestamps,np.zeros_like(timestamps)+count, c=color, marker=marker , s=2)
        
        if j%2==0:
            for spike in timestamps:
                if spike in prev_timestamps:
                    plt.vlines(spike, count-20, count+20, color='k', linewidth=4)
        
        prev_timestamps=timestamps
    plt.xlim(0,1000)
    plt.ylim(0,count)
    plt.xticks([0,250,500,750,999], ['-1','-0.75','-0.5','-0.25','0'])
    plt.yticks(np.arange(0,count,300), ['0','50','100','150','200','250', '300','350','400','450'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time (secons)')
    plt.ylabel('Trials')
    
    #Zoom
    fig,ax= plt.subplots(1,1, figsize=(4,8))  
    count=0
    prev_timestamps=[0,0]
    for j,row in enumerate(Spikes_NeuronPair):
        if j%2==0:
            color='c'
            marker='.'
            count+=2
        else:
            color='m'
            count+=4
            marker='.'
        timestamps=[i for i,x in enumerate(row) if x]
        plt.scatter(timestamps,np.zeros_like(timestamps)+count, c=color, marker=marker , s=10)
        
        if j%2==0:
            for spike in timestamps:
                if spike in prev_timestamps:
                    plt.vlines(spike, count-2, count, color='k', linewidth=4)
        
        prev_timestamps=timestamps
    plt.xlim(0,1000)
    plt.ylim(210*6,(240*6)+1)
    plt.xticks([0,250,500,750,999], ['-1','-0.75','-0.5','-0.25','0'])
    plt.yticks(np.arange(1260,1441,60), ['210','220','230','240'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time (secons)')
    plt.ylabel('Trials')
    return

###############################################################################
# XCorr_example
###############################################################################
    
def XCorr_example(XCorr_df, pair):#uncorrelated: 10, correlated:23
    Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
    Filter=np.array([0.05,0.25,0.4,0.25,0.05])
    Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter)
    plt.figure()
    plt.plot(np.arange(1004),Subtracted_XCorr, color='k',alpha=0.5)
    plt.xticks([0,251,502, 753, 1004], ['-1', '-0.5','0','0.5','1'])
    plt.yticks( [0.02,0.04,0.06,0.08], [x*2.5 for x in [0.02,0.04,0.06,0.08]])
    plt.title(XCorr_df.loc[pair,'Neuron1_name'] +'\n'+ XCorr_df.loc[pair,'Neuron2_name'])
    
    #zoom inset: no filter
    Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
    plt.figure()
    plt.plot(np.arange(31),Subtracted_XCorr[484:515], color='k',alpha=0.5) #the exact 0-lag is index 499
    plt.xticks([5,15,25], ['-10','0','10'])
    plt.yticks( [0,0.04,0.08, 0.12])
    plt.vlines(15, 0,0.12, color='grey', linestyles='dotted')
    plt.title(XCorr_df.loc[pair,'Neuron1_name'] +'\n'+ XCorr_df.loc[pair,'Neuron2_name'] + ' NO FILTER')
    
    return

###############################################################################
# Claustrum_v_S1_XCorr
###############################################################################
    
def Claustrum_v_S1_XCorr(XCorr_df, S1XCorr_df):
    from scipy.stats import chi2_contingency 
    
    Corr_LickPrefs1=[]
    Corr_LickPrefs2=[]
    Corr_strength=[]
    FR=[]
    FRtype=[]
    Test_group=[]
    Neuron1_name=[]
    Neuron2_name=[]
    Neuron1_Category=[]
    Neuron2_Category=[]
    optocount=0
    count=0
    #Sig_index=[]
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
            FR.append(XCorr_df.loc[pair,'Geometric mean FR'])
            FRtype.append(XCorr_df.loc[pair,'Geometric mean FR LEFT'])
            Neuron1_name.append(XCorr_df.loc[pair,'Neuron1_name'])
            Neuron2_name.append(XCorr_df.loc[pair,'Neuron2_name'])
            Neuron1_Category.append(XCorr_df.loc[pair,'Neuron1_Category'])
            Neuron2_Category.append(XCorr_df.loc[pair,'Neuron2_Category'])
            if ((XCorr_df.loc[pair,'Neuron1_Category'] == 'OptoTag' ) & (XCorr_df.loc[pair,'Neuron2_Category']=='OptoTag')):
                optocount+=1
                #print(count)
        count+=1

    S1Corr_LickPrefs1=[]
    S1Corr_LickPrefs2=[]
    S1Corr_strength=[]
    S1FR=[]
    S1FRtype=[]
    S1mouse=[]
    #S1Sig_index=[]
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
            S1FR.append(S1XCorr_df.loc[pair,'Geometric mean FR'])
            S1FRtype.append(S1XCorr_df.loc[pair,'Geometric mean FR LEFT'])
            S1mouse.append(S1XCorr_df.loc[pair,'Neuron1_name'])

    Sig=Corr_strength
    FRSig=FR
    LickPref1_Sig= Corr_LickPrefs1
    LickPref2_Sig= Corr_LickPrefs2
    S1Sig=S1Corr_strength
    S1FRSig=S1FR
    
    All_count=0
    for pair in XCorr_df.index:
        if ( (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name'])):
            All_count+=1
    S1All_count=0
    for pair in S1XCorr_df.index:
        if ( (S1XCorr_df.loc[pair,'Neuron1_name'] != S1XCorr_df.loc[pair,'Neuron2_name'])):
            S1All_count+=1
            
    plt.figure()
    plt.pie([len(Sig), All_count-len(Sig)], 
            colors=['k', 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[len(Sig), All_count-len(Sig)],
            labeldistance=0.65 )
    plt.title('Claustrum')        
    
    plt.figure()
    plt.pie([len(S1Sig), S1All_count-len(S1Sig)], 
            colors=['k', 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[len(S1Sig), S1All_count-len(S1Sig)],
            labeldistance=0.65 )
    plt.title('S1') 
    
    g, p, dof, expctd=chi2_contingency([[len(Sig),len(S1Sig)],[All_count-len(Sig),S1All_count-len(S1Sig)]])

        
    return p, Sig, S1Sig, LickPref1_Sig, LickPref2_Sig, FRSig, S1FRSig

###############################################################################
# Claustrum_v_S1_XCorr
###############################################################################
    
def CumFreq_Claustrum_v_S1(Sig, S1Sig):
    plt.style.use('default')
    RR=Sig
    RR=[np.log(x*100/0.4) for x in RR]
    Cum=stats.cumfreq(RR, numbins=100)#, defaultreallimits=(-0.5,1))
    x=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    RRplot=plt.plot(x,Cum.cumcount/np.size(RR), color='k', linestyle='dotted')
    LL=S1Sig 
    LL=[np.log(x*100/0.4) for x in LL]
    Cum=stats.cumfreq(LL, numbins=100)#, defaultreallimits=(-0.5,1))
    x=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    LLplot=plt.plot(x,Cum.cumcount/np.size(LL), color='grey', linestyle='dotted')
    plt.xticks([np.log(x) for x in [1,2,4,8,16,32, 64]], [1,2,4,8,16,32, 64])
    plt.xlabel('Percentage coincident spikes')
    plt.ylabel('ECDF')
    plt.legend(['Claustrum, n='+str(len(Sig)),'S1, n='+str(len(S1Sig))])
    s,p=stats.mannwhitneyu(Sig, S1Sig)
    plt.title('MannWhitneyU p='+str(p))
    
    Sig_mean=np.mean(Sig)
    Sig_sem=np.std(Sig)/np.sqrt(len(Sig))
    S1Sig_mean=np.mean(S1Sig)
    S1Sig_sem=np.std(S1Sig)/np.sqrt(len(S1Sig))
    return Sig_mean, Sig_sem, S1Sig_mean, S1Sig_sem

###############################################################################
# Quandrant_plot
###############################################################################
    
def Quandrant_plot(XCorr_df, Sig, LickPref1_Sig, LickPref2_Sig, cutoff):
    
    Corr_LickPrefs1=[]
    Corr_LickPrefs2=[]
    Corr_strength=[]
    FR=[]
    FRtype=[]
    Test_group=[]
    Neuron1_name=[]
    Neuron2_name=[]
    #Sig_index=[]
    for pair in XCorr_df.index:
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
        Filter=np.array([0.05,0.25,0.4,0.25,0.05])
        Subtracted_XCorr=np.convolve(Subtracted_XCorr,Filter) 
        Baseline=( np.mean(Subtracted_XCorr[:100]) + np.mean(Subtracted_XCorr[-100:]) ) / 2
        if ( (XCorr_df.loc[pair,'Neuron1_name'] != XCorr_df.loc[pair,'Neuron2_name'])):
            Corr_LickPrefs1.append(XCorr_df.loc[pair,'Neuron1_LickPref'])
            Corr_LickPrefs2.append(XCorr_df.loc[pair,'Neuron2_LickPref'])
            Corr_strength.append( np.max(Subtracted_XCorr[488:512]) - Baseline )  
            FR.append(XCorr_df.loc[pair,'Geometric mean FR'])
            FRtype.append(XCorr_df.loc[pair,'Geometric mean FR LEFT'])
            Neuron1_name.append(XCorr_df.loc[pair,'Neuron1_name'])
            Neuron2_name.append(XCorr_df.loc[pair,'Neuron2_name'])
    All=Corr_strength
    FRAll=FR
    LickPref1_All= Corr_LickPrefs1
    LickPref2_All= Corr_LickPrefs2
    
    plt.style.use('default')
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    plt.fill([-1,-0.1,-0.1,-1], [-1,-1,-0.1,-0.1], color='#00D88A', alpha=0.1)
    plt.fill([0.1,1,1,0.1], [0.1,0.1,1,1], color='#A1007D', alpha=0.1)
    plt.fill([-1,-0.1,-0.1,-1], [0.1,0.1,1,1], color='grey', alpha=0.1)
    plt.fill([0.1,1,1,0.1], [-1,-1,-0.1,-0.1], color='grey', alpha=0.1)
    
    NoOutlyer=[x if x>0.0001 else 0.0001 for x in All]    
    NoOutlyer=[np.log(x) for x in NoOutlyer]
    plt.scatter(LickPref1_All,LickPref2_All, c='grey', alpha=0.1, linewidth=0)
    NoOutlyer=[x if x>0.0001 else 0.0001 for x in Sig]    
    NoOutlyer=[np.log(x) for x in NoOutlyer]
    plt.scatter(LickPref1_Sig,LickPref2_Sig, c='k', linewidth=0)
    plt.scatter([y for x,y,z in zip(Sig, LickPref1_Sig, LickPref2_Sig) if ( ((y>cutoff) & (z>cutoff)) | ((y<-cutoff) & (z<-cutoff)) | ((y>cutoff) & (z<-cutoff)) | ((y<-cutoff) & (z>cutoff)) )],
                 [z for x,y,z in zip(Sig, LickPref1_Sig, LickPref2_Sig) if ( ((y>cutoff) & (z>cutoff)) | ((y<-cutoff) & (z<-cutoff)) | ((y>cutoff) & (z<-cutoff)) | ((y<-cutoff) & (z>cutoff)) )])
    plt.vlines(0,-1,1)
    plt.hlines(0,-1,1)
    
    return All, FRAll, LickPref1_All, LickPref2_All


###############################################################################
# PiePlot_LickPref
###############################################################################
def PiePlot_LickPref(Sig,LickPref1_Sig, LickPref2_Sig, LickPref1_All, LickPref2_All):
    from scipy.stats import chi2_contingency 
    cutoff=0.1
    RR=sum([1 for x,y in zip(LickPref1_All,LickPref2_All) if ((x>cutoff) & (y>cutoff))])
    LL=sum([1 for x,y in zip(LickPref1_All,LickPref2_All) if ((x<-cutoff) & (y<-cutoff))])
    RL=sum([1 for x,y in zip(LickPref1_All,LickPref2_All) if ( ((x<-cutoff) & (y>cutoff)) | ((x>cutoff) & (y<-cutoff)) )])
    plt.figure()
    plt.pie([RR,LL,RL], 
            colors=['#A1007D', '#00D88A', 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[RR,LL,RL],
            labeldistance=0.65 )
    plt.title('Composition of pairs(0.1 cutoff)')       
    
    RRSig=len([x for x,y,z in zip(Sig, LickPref1_Sig, LickPref2_Sig) if ( (y>cutoff) & (z>cutoff))]  ) 
    LLSig=len([x for x,y,z in zip(Sig, LickPref1_Sig, LickPref2_Sig) if ( (y<-cutoff) & (z<-cutoff))]  )
    RLSig=len([x for x,y,z in zip(Sig, LickPref1_Sig, LickPref2_Sig) if ( ( (y>cutoff) & (z<-cutoff)) |  ( (y<-cutoff) & (z>cutoff)) )]  )
    plt.figure()
    plt.pie([RRSig,LLSig,RLSig], 
            colors=['#A1007D', '#00D88A', 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[RRSig,LLSig,RLSig],
            labeldistance=0.65 )
    plt.title('Composition of SIG pairs(0.1 cutoff)')  
    
    g, p, dof, expctd=chi2_contingency([[RR-RRSig,LL-LLSig,RL-RLSig],[RRSig, LLSig, RLSig]])
    
    return p