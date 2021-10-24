###############################################################################
# Figure S6
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
    # - Compare_XCorrWidth_Claustrum_vs_S1()
    # - Claustrum_v_S1_XCorr_poststim()

###############################################################################
# Figure S6C-H
###############################################################################
    
def Compare_XCorrWidth_Claustrum_vs_S1(XCorr_df, S1XCorr_df):
    # examples for CLaustrum 1878, 1235
    Corr_LickPrefs1=[]
    Corr_LickPrefs2=[]
    Corr_strength=[]
    FR=[]
    FRtype=[]
    Neuron1_name=[]
    Neuron2_name=[]
    Neuron1_Category=[]
    Neuron2_Category=[]
    
    peak_indices=[]
    widths=[]
    widths2=[]
    count=0
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
            Corr_LickPrefs1.append(XCorr_df.loc[pair,'Neuron1_LickPref'])
            Corr_LickPrefs2.append(XCorr_df.loc[pair,'Neuron2_LickPref'])
            Corr_strength.append( np.max(Subtracted_XCorr[488:512]) - Baseline )  
            FR.append(XCorr_df.loc[pair,'Geometric mean FR'])
            FRtype.append(XCorr_df.loc[pair,'Geometric mean FR LEFT'])
            Neuron1_name.append(XCorr_df.loc[pair,'Neuron1_name'])
            Neuron2_name.append(XCorr_df.loc[pair,'Neuron2_name'])
            Neuron1_Category.append(XCorr_df.loc[pair,'Neuron1_Category'])
            Neuron2_Category.append(XCorr_df.loc[pair,'Neuron2_Category'])
            
            #Find width of peak
            #get peak, get heights, find indices of halfmax
            Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
            peak_idx=np.argmax(Subtracted_XCorr)
#            peak_val=Subtracted_XCorr[peak_idx]-Baseline
#            half_peak_val=Baseline+ (peak_val/2)
#            start_idx=peak_idx-[i for i,x in enumerate(Subtracted_XCorr[peak_idx:0:-1]) if ((x<half_peak_val) & (Subtracted_XCorr[peak_idx-i-1]<half_peak_val) & (Subtracted_XCorr[peak_idx-i-2]<half_peak_val))][0]
#            stop_idx=peak_idx+[i for i,x in enumerate(Subtracted_XCorr[peak_idx:-3]) if ((x<half_peak_val) & (Subtracted_XCorr[peak_idx+i+1]<half_peak_val) & (Subtracted_XCorr[peak_idx+i+2]<half_peak_val))][0]
#            peak_indices.append(peak_idx)
#            widths.append(stop_idx-start_idx-1)
            
             # Alternative strategy:
            Rise_onset=peak_idx-[i for i,x in enumerate(Subtracted_XCorr[peak_idx:0:-1]) if ((x<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx-i-1]<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx-i-2]<(Baseline+3*BaselineSD)))][0]
            Fall_offset=peak_idx+[i for i,x in enumerate(Subtracted_XCorr[peak_idx:-3]) if ((x<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx+i+1]<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx+i+2]<(Baseline+3*BaselineSD)))][0]
            widths2.append(Fall_offset-Rise_onset-1)
            if count in [1878, 1235]:                                        
                Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
                plt.figure()
                plt.plot(np.arange(200),Subtracted_XCorr[400:600], color='k',alpha=0.5) #the exact 0-lag is index 502
                plt.xticks([0,24,49,74,99,124,149,174,199], ['-100','-75','-50','-25','0','25','50','75','100'])
                plt.vlines(peak_idx-400, 0,0.02, color='grey', linestyles='dotted')
                plt.vlines(Fall_offset-400, 0,0.02, color='grey', linestyles='dotted')
                plt.vlines(Rise_onset-400, 0,0.02, color='grey', linestyles='dotted')
                plt.title(XCorr_df.loc[pair,'Neuron1_name'] +'\n'+ XCorr_df.loc[pair,'Neuron2_name'] + ' \nWidth=' + str(Fall_offset-Rise_onset) + ' ms\nCount='+str(count))
                plt.xlabel('Lag (ms)')
                plt.ylabel('Fraction of coincident spikes')
                plt.ylim(-0.009,0.03)
        count+=1
    Fig,ax=plt.subplots(1,1, figsize=(3,5))
    plt.hist([x-499 for x in peak_indices],bins=7, range=[-3.5,3.5], align='mid')
    plt.xlabel('Lag (ms)')
    plt.ylabel('Counts')
    plt.title('Claustrum')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    #S1 examples: 30, 100, 107
    S1Corr_LickPrefs1=[]
    S1Corr_LickPrefs2=[]
    S1Corr_strength=[]
    S1FR=[]
    S1FRtype=[]
    
    S1peak_indices=[]
    S1widths=[]
    S1widths2=[]
    count=0
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
            S1Corr_LickPrefs1.append(S1XCorr_df.loc[pair,'Neuron1_LickPref'])
            S1Corr_LickPrefs2.append(S1XCorr_df.loc[pair,'Neuron2_LickPref'])
            S1Corr_strength.append( np.max(Subtracted_XCorr[488:512]) - Baseline )       
            S1FR.append(S1XCorr_df.loc[pair,'Geometric mean FR'])
            S1FRtype.append(S1XCorr_df.loc[pair,'Geometric mean FR LEFT'])
            
            #Find width of peak
            #get peak, get heights, find indices of halfmax
            Subtracted_XCorr=S1XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-S1XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
            peak_idx=np.argmax(Subtracted_XCorr)
#            peak_val=Subtracted_XCorr[peak_idx]-Baseline
#            half_peak_val=Baseline+ (peak_val/2)
#            start_idx=peak_idx-[i for i,x in enumerate(Subtracted_XCorr[peak_idx:0:-1]) if ((x<half_peak_val) & (Subtracted_XCorr[peak_idx-i-1]<half_peak_val) & (Subtracted_XCorr[peak_idx-i-2]<half_peak_val))][0]
#            stop_idx=peak_idx+[i for i,x in enumerate(Subtracted_XCorr[peak_idx:-3]) if ((x<half_peak_val) & (Subtracted_XCorr[peak_idx+i+1]<half_peak_val) & (Subtracted_XCorr[peak_idx+i+2]<half_peak_val))][0]
#            S1peak_indices.append(peak_idx)
#            S1widths.append(stop_idx-start_idx-1)
            
            # Alternative strategy:
            S1Rise_onset=peak_idx-[i for i,x in enumerate(Subtracted_XCorr[peak_idx:0:-1]) if ((x<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx-i-1]<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx-i-2]<(Baseline+3*BaselineSD)))][0]
            S1Fall_offset=peak_idx+[i for i,x in enumerate(Subtracted_XCorr[peak_idx:-3]) if ((x<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx+i+1]<(Baseline+3*BaselineSD)) & (Subtracted_XCorr[peak_idx+i+2]<(Baseline+3*BaselineSD)))][0]
            S1widths2.append(S1Fall_offset-S1Rise_onset-1)
            if count in [30,100,107]:
 
                Subtracted_XCorr=S1XCorr_df.loc[pair,'XCorr-50to50_ALL1sec']-S1XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1sec']
                plt.figure()
                plt.plot(np.arange(200),Subtracted_XCorr[400:600], color='k',alpha=0.5) #the exact 0-lag is index 502
                plt.xticks([0,24,49,74,99,124,149,174,199], ['-100','-75','-50','-25','0','25','50','75','100'])
                plt.vlines(peak_idx-400, 0,0.02, color='grey', linestyles='dotted')
                plt.vlines(S1Fall_offset-400, 0,0.02, color='grey', linestyles='dotted')
                plt.vlines(S1Rise_onset-400, 0,0.02, color='grey', linestyles='dotted')
                plt.title(S1XCorr_df.loc[pair,'Neuron1_name'] +'\n'+ S1XCorr_df.loc[pair,'Neuron2_name'] + ' \nWidth=' + str(S1Fall_offset-S1Rise_onset) + ' ms\nCount='+str(count))
                plt.xlabel('Lag (ms)')
                plt.ylabel('Fraction of coincident spikes')
                plt.ylim(-0.009,0.03)        
            count+=1
        
#    Fig,ax=plt.subplots(1,1, figsize=(4,5))
#    plt.hist([x-499 for x in peak_indices],bins=21, range=[-10.5,10.5], align='mid', alpha=0.5)
#    plt.hist([x-499 for x in S1peak_indices],bins=21, range=[-10.5,10.5], align='mid', alpha=0.5)
#    plt.legend(['Claustrum','S1'])
#    plt.xlabel('Lag (ms)')
#    plt.ylabel('Counts')
#    plt.title('Distribution of CCG Lag')
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)        
#    
#    Fig,ax=plt.subplots(1,1, figsize=(4,5))
#    plt.hist([x-499 for x in peak_indices],bins=21, range=[-10.5,10.5], align='mid', alpha=0.5, stacked=True, density=True)
#    plt.hist([x-499 for x in S1peak_indices],bins=21, range=[-10.5,10.5], align='mid', alpha=0.5, stacked=True, density=True)
#    plt.legend(['Claustrum','S1'])
#    plt.xlabel('Lag (ms)')
#    plt.ylabel('Proportion')
#    plt.ylim(0,1)
#    s,p=stats.mannwhitneyu(peak_indices, S1peak_indices)
#    plt.title('p='+str(p))
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)      
#    
#    Fig,ax=plt.subplots(1,1, figsize=(4,5))
#    plt.hist([x for x in widths2],bins=50, range=[0.5,80.5], align='mid', alpha=0.5)
#    plt.hist([x for x in S1widths2],bins=50, range=[0.5,80.5], align='mid', alpha=0.5)
#    plt.legend(['Claustrum','S1'])
#    plt.xlabel('width (ms)')
#    plt.ylabel('Counts')
#    plt.title('Distribution of CCG width')
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)       
    
    Fig,ax=plt.subplots(1,1, figsize=(4,5))
    plt.hist([x for x in widths2],bins=50, range=[0.5,80.5], align='mid', alpha=0.5, stacked=True, density=True)
    plt.hist([x for x in S1widths2],bins=50, range=[0.5,80.5], align='mid', alpha=0.5, stacked=True, density=True)
    plt.legend(['Claustrum','S1'])
    plt.xlabel('width (ms)')
    plt.ylabel('Proportion')
    plt.ylim(0,0.4)
    s,p=stats.mannwhitneyu(widths2, S1widths2)
    plt.title('p='+str(p))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)      
    
    mean_w=np.mean(widths2)
    sem_w=np.std(widths2)/np.sqrt(len(widths2))
    mean_S1w=np.mean(S1widths2)
    sem_S1w=np.std(S1widths2)/np.sqrt(len(S1widths2))
    #Results:mean_w:3.493 , sem_w:0.486 , mean_S1w:21.597 , sem_S1w:2.210, n_w=75 , n_S1w=129
    
    Sig=Corr_strength
    FRSig=FR
    LickPref1_Sig= Corr_LickPrefs1
    LickPref2_Sig= Corr_LickPrefs2
    S1Sig=S1Corr_strength
    S1FRSig=S1FR
    return mean_w, sem_w, mean_S1w, sem_S1w, p



###############################################################################
# Figure S6A,B
###############################################################################
def Claustrum_v_S1_XCorr_poststim(XCorr_df, S1XCorr_df):
    from scipy.stats import chi2_contingency 

    Corr_LickPrefs1=[]
    Corr_LickPrefs2=[]
    Corr_strength=[]
    FR=[]
    Neuron1_name=[]
    Neuron2_name=[]
    Neuron1_Category=[]
    Neuron2_Category=[]
    count=0
    for pair in XCorr_df.index:
        
        Subtracted_XCorr=XCorr_df.loc[pair,'XCorr-50to50_ALL1secPOST']-XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1secPOST']
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
            FR.append(XCorr_df.loc[pair,'Geometric mean FR ALLPOST'])
            Neuron1_name.append(XCorr_df.loc[pair,'Neuron1_name'])
            Neuron2_name.append(XCorr_df.loc[pair,'Neuron2_name'])
            Neuron1_Category.append(XCorr_df.loc[pair,'Neuron1_Category'])
            Neuron2_Category.append(XCorr_df.loc[pair,'Neuron2_Category'])
#            plt.figure()
#            plt.plot(np.arange(200),Subtracted_XCorr[400:600], color='k',alpha=0.5) #the exact 0-lag is index 502
                #print(count)
        count+=1
    
    S1Corr_LickPrefs1=[]
    S1Corr_LickPrefs2=[]
    S1Corr_strength=[]
    S1FR=[]
    for pair in S1XCorr_df.index:
        Subtracted_XCorr=S1XCorr_df.loc[pair,'XCorr-50to50_ALL1secPOST']-S1XCorr_df.loc[pair,'XCorr-50to50_shuffle_ALL1secPOST']
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
            S1FR.append(S1XCorr_df.loc[pair,'Geometric mean FR ALLPOST'])
#            plt.figure()
#            plt.plot(np.arange(200),Subtracted_XCorr[400:600], color='k',alpha=0.5) #the exact 0-lag is index 502
                #print(count)
    
    Sig=Corr_strength
    S1Sig=S1Corr_strength
    
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
    plt.title('Claustrum POST1sec')        
    
    plt.figure()
    plt.pie([len(S1Sig), S1All_count-len(S1Sig)], 
            colors=['k', 'darkgrey'], 
            startangle=90 , 
            counterclock=False, 
            labels=[len(S1Sig), S1All_count-len(S1Sig)],
            labeldistance=0.65 )
    plt.title('S1 POST1sec, chisquare p=0.0193') 
    
    g, p, dof, expctd=chi2_contingency([[len(Sig),len(S1Sig)],[All_count-len(Sig),S1All_count-len(S1Sig)]])
    
    plt.figure()
    plt.style.use('default')
    RR=Sig#[x for x in Sig if (x>0.0001)]  
    RR=[np.log(x*100/0.4) for x in RR]
    Cum=stats.cumfreq(RR, numbins=100)#, defaultreallimits=(-0.5,1))
    x=Cum.lowerlimit + np.linspace(0, Cum.binsize*Cum.cumcount.size, Cum.cumcount.size)
    RRplot=plt.plot(x,Cum.cumcount/np.size(RR), color='k', linestyle='dotted')
    LL=S1Sig #[x for x in S1Sig if (x>0.0001)]    
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
    
    
    Sig_mean=np.mean(Sig)
    Sig_sem=np.std(Sig)/np.sqrt(len(Sig))
    S1Sig_mean=np.mean(S1Sig)
    S1Sig_sem=np.std(S1Sig)/np.sqrt(len(S1Sig))
    #REsults: Sig_mean=0.0678 , Sig_sem=0.0111 , S1Sig_mean=0.0550 , S1Sig_sem=0.00761
    return Sig_mean, Sig_sem, S1Sig_mean, S1Sig_sem, p
