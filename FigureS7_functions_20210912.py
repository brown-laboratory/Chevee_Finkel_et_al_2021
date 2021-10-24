###############################################################################
# Figure S7
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
    # - clustering()
    # - Clusters_bargraph()

###############################################################################
# Figure S7A,B
###############################################################################

def clustering(mice, master_log):
    from scipy.stats.stats import pearsonr
    import itertools
    from numpy import ma
    from matplotlib import cbook
    from matplotlib import colors
    from matplotlib.colors import Normalize
    
    class MidpointNormalize(colors.Normalize):
    	"""
    	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
    
    	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    	"""
    	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
    		self.midpoint = midpoint
    		colors.Normalize.__init__(self, vmin, vmax, clip)
    
    	def __call__(self, value, clip=None):
    		# I'm ignoring masked values and all kinds of edge cases to make a
    		# simple example...
    		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
    		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    # for each neuron, make an 1d array with all 8 trial type in a row. Seconds 0 to 1 in 25msec bins
    Correlations=pd.DataFrame() #initiate a df
    Name_df=['unit_name','OnedVector', 'Category']
    for i,name in enumerate(Name_df):   
        Correlations[Name_df[i]]=np.zeros(len(Correlations))
        Correlations[Name_df[i]]=Correlations[Name_df[i]].astype(object)
        
    
    ## Build input matrix to pairwise correlation
    #OPTION 1: Align by modality
    counter=0
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                print(unit)
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                
                OnedAllTrials=np.array(0) #initialize a list to put all trial together
                for trialtype in [ 'SomHit', 'VisHit']:#, 'SomCR', 'VisCR']:  #IMPORTANT: CHOICE OF WHAT TRIALS TO CUSTER ON
                    unit_log_trialtype=unit_log[unit_log['Stim/Block/Response']==trialtype] #get only the trials from the proper type
                    if not unit_log_trialtype.empty:
                        Onesec_25msecbins_mean=np.mean(unit_log_trialtype['Zscored_-1to3sec_25msecbins_StimAligned'])[:159]  #get the one sec ([40:80]) follwing stim 
                        OnedAllTrials=np.hstack((OnedAllTrials, Onesec_25msecbins_mean)) #stick it to the 1d vector  
#                
                Correlations.at[counter, 'unit_name']=[mouse + '_' + day[0] + '_' + unit[0]]
                Correlations.at[counter, 'OnedVector']=OnedAllTrials
                Correlations.at[counter, 'Category']=np.unique(unit_log['Category'])[0]
                counter+=1
   
    # repeat the same loop but store the results in a 2d matrix to then plot a similarity matrix
    correlation_matrix=np.zeros((len(Correlations), len(Correlations)))
    count1=0;#have to use counters because the indeces are off dur to the dropped sessions in paragraph above ## DEPRECATED (but doesn't interfere with code)
    for each1 in Correlations.index:
        count2=0
        for each2 in Correlations.index:
            correlation_matrix[count1,count2] = pearsonr(Correlations.loc[each1,'OnedVector' ],Correlations.loc[each2,'OnedVector' ])[0] # add [0] to only get the r amd not the p-value
            count2+=1
        count1+=1
        
    df=pd.DataFrame(data=correlation_matrix, columns=Correlations['unit_name'])#create df to input clustermap
    
    lut = dict(zip(Correlations['Category'].unique(), "rbg"))
    Correlations_Category_indexReset=  Correlations['Category'].reset_index() #need to reset the index because otherwise it's off because of the nans we dropped earlier ## DEPRECATED (but doesn't interfere with code)
    row_colors =Correlations_Category_indexReset['Category'].map(lut) #a list of colors coresponding to 'Category'
    
    ## Cluster and draw clustermap
    from scipy.spatial import distance
    from scipy.cluster import hierarchy
    Z=hierarchy.linkage(correlation_matrix,method='average') #do it once to get Z which you will use below to get clusters
    row_linkage = hierarchy.linkage( #these you're only doing for the plot
        distance.pdist(correlation_matrix), method='average')
    col_linkage = hierarchy.linkage( #these you're only doing for the plot
        distance.pdist(correlation_matrix.T), method='average')
    
#    sns.set(font_scale=1)
#    map_test=sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage,  method="average", xticklabels=1, figsize=(15,15), annot_kws={"size": 6}, row_colors=row_colors)
#    map_test.ax_col_dendrogram.set_xlim([0,0]) #removes top dendogram on clustemap
    
#    # Draw the legend bar for the classes       (have to do a trick here: draw an invisible bar plot and then plot the legened for that)          
#    for label in Correlations['Category'].unique():
#        map_test.ax_col_dendrogram.bar(0, 0, color=lut[label],label=label, linewidth=0)
#    map_test.ax_col_dendrogram.legend(loc="center", ncol=5)
    
    ## Draw dendogram
    plt.figure()
    Cluster_assignment=hierarchy.fcluster(Z, 5, criterion='distance')
    D=hierarchy.dendrogram(Z, color_threshold=5) #13 groups
    #plt.savefig('Dendogram_6groups_5p2_ALL.pdf')
    

    
    cluster_colors=['k','c','m','darkorange','r','g','b','lime', 'gold', 'purple','olivedrab', 'k', 'c', 'm']
    #Option #1
    temp_units1=[]
    mouse_name=[]
    Opto=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                mouse_name.append(mouse)
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                if unit_log['Category'].values[0]=='OptoTag':
                    Opto.append(1)
                else:
                        Opto.append(0)
                norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
    #            if mouse in ORIGINAL:
                temp_trial=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    #            elif mouse in REVERSED:
    #                temp_trial=[]
    #                for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
    #                    trial_log=unit_log.loc[trial,:]
    #                    temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
                temp_units1.append(np.mean(temp_trial, axis=0)/norm)
    
    temp_units2=[]
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                norm=np.mean(np.mean(unit_log['-1to3sec_25msecbins_StimAligned'],0)[:39])
    #            if mouse in ORIGINAL:
                temp_trial=[]
                for trial in unit_log[unit_log['Stim/Block/Response']=='VisHit'].index:
                    trial_log=unit_log.loc[trial,:]
                    temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
    #            elif mouse in REVERSED:
    #                temp_trial=[]
    #                for trial in unit_log[unit_log['Stim/Block/Response']=='SomHit'].index:
    #                    trial_log=unit_log.loc[trial,:]
    #                    temp_trial.append(trial_log['-1to3sec_25msecbins_StimAligned']) 
                temp_units2.append(np.mean(temp_trial, axis=0)/norm)
    
                
    sorted_temp_units1=[temp_units1[row] for row in D['leaves']]
    sorted_temp_units2=[temp_units2[row] for row in D['leaves']]
    sorted_mouse_name=[mouse_name[row] for row in D['leaves']]
    sorted_opto=[Opto[row] for row in D['leaves']]
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plt.style.use('default')
    fig, ax=plt.subplots(1,2,figsize=(6,8))    
    
    SomHit=ax[0].imshow(sorted_temp_units1,norm=MidpointNormalize(midpoint=1.,vmin=0, vmax=4), cmap='bwr')
    ax[0].get_yaxis().set_visible(False)
    
    #colorbar with cluster
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    sorted_colors=[row_colors[each] for each in D['leaves']] #'Category' colors
    sorted_cluster=[Cluster_assignment[each] for each in D['leaves']] # Cluster_group
    for i,each in enumerate(reversed(sorted_colors)):
        cax.bar(1,1,bottom=i, color=cluster_colors[sorted_cluster[np.size(sorted_temp_units1,0)-1-i]], edgecolor='none') # reverse the sort because we are plotting from top to bottom
    cax.set_ylim(0,np.size(sorted_temp_units1,0))
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    
    #colorbar with ORIGINAL vs REVERSED
    cax = divider.append_axes("left", size="5%", pad=0.05)
    for i,mouse in enumerate(reversed(sorted_mouse_name)):
        if mouse in ORIGINAL:
            cax.bar(1,1,bottom=i, color='k', edgecolor='none') # reverse the sort because we are plotting from top to bottom
        elif mouse in REVERSED:
            cax.bar(1,1,bottom=i, color='grey', edgecolor='none') # reverse the sort because we are plotting from top to bottom
    cax.set_ylim(0,np.size(sorted_temp_units1,0))
    cax.get_xaxis().set_visible(False)
    
    
    
    #horizontal lines
    h_pos=[i  for i,x in enumerate(sorted_cluster[1:]) if x>sorted_cluster[i]]
    ax[0].hlines(h_pos,0,158, color='k', linewidth=2)
    ax[0].vlines(39.5,0,len(temp_units1)-1, color='k', linewidth=2)
    ax[0].set_title('SomHit')
    ax[0].set_xticks([-0.5, 39.5,79.5,119.5, 159.5])
    ax[0].set_xticklabels( [ '-1','0', '1', '2', '3'])
    
    SomHit=ax[1].imshow(sorted_temp_units2,norm=MidpointNormalize(midpoint=1.,vmin=0, vmax=4), cmap='bwr')
    
    #color scale
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(SomHit, cax=cax)
    
    #horizontal lines
    sorted_colors=[row_colors[each] for each in D['leaves']] #'Category' colors
    sorted_cluster=[Cluster_assignment[each] for each in D['leaves']] # Cluster_group
    h_pos=[i  for i,x in enumerate(sorted_cluster[1:]) if x>sorted_cluster[i]]
    ax[1].hlines(h_pos,0,158, color='k', linewidth=2)
    ax[1].vlines(39.5,0,len(temp_units2)-1, color='k', linewidth=2)
    ax[1].set_title('VisHit')
    ax[1].set_xticks([-0.5, 39.5,79.5,119.5, 159.5])
    ax[1].set_xticklabels( [ '-1','0', '1', '2', '3'])
    
    #colorbar with OptoTag
    opto_colors=['yellow','purple']
    cax = divider.append_axes("left", size="5%", pad=0.05)
    sorted_colors=[opto_colors[each] for each in sorted_opto] #'Category' colors
    for i,each in enumerate(reversed(sorted_colors)):
        cax.bar(1,1,bottom=i, color=opto_colors[sorted_opto[np.size(sorted_temp_units1,0)-1-i]], edgecolor='none') # reverse the sort because we are plotting from top to bottom
    cax.set_ylim(0,np.size(sorted_temp_units1,0))
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    
    fig.suptitle('ALL')
    
    
    
    ##########
    # AUC MEANS
    ##########
    Cluster_Category=[1,2,3,4,5,6,7,8,9,10, 11, 12, 13]
    #Cluster_Category=[1,2,3]
    cluster_colors=['k','c','m','darkorange','r','g','b','lime', 'gold', 'purple','olivedrab', 'k', 'k','c']
    
    plt.style.use('default')
    fig, ax=plt.subplots(13,1,figsize=(4,25))
    #fig, ax=plt.subplots(3,1,figsize=(5,10))
    for i in Cluster_Category:
        temp_units1=[]
        temp_units2=[]
        Touchcount=0
        Visualcount=0
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name'][(day_log['test']==i)]):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    if mouse in ORIGINAL:
                        temp_units1.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
                        temp_units2.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
                    if mouse in REVERSED:
                        temp_units1.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisCR'].values[0])
                        temp_units2.append(unit_log['unit_mean_25msecAUC_-1to3sec_stimAligned_VisHit_vs_SomCR'].values[0])
       
                    
       
    #    if not [x for x in temp_units1]:
    #        temp_units1=np.zeros((2,159))
    #    if not [x for x in temp_units2]:
    #        temp_units2=np.zeros((2,159))
    #  
        #Calculate/plot Mean/SEM
        SEM1=sp.stats.sem(temp_units1, axis=0)
        SomHit, =ax[i-1].plot(np.arange(159), np.mean(temp_units1, axis=0), color=cluster_colors[i])
        ax[i-1].fill_between(np.arange(159),np.mean(temp_units1, axis=0)+SEM1, np.mean(temp_units1, axis=0)-SEM1, color=cluster_colors[i], alpha=0.1)
        SEM2=sp.stats.sem(temp_units2, axis=0)
        VisHit, =ax[i-1].plot(np.arange(159), np.mean(temp_units2, axis=0), color='darkgrey')
        ax[i-1].fill_between(np.arange(159),np.mean(temp_units2, axis=0)+SEM2, np.mean(temp_units2, axis=0)-SEM2, color='darkgrey', alpha=0.1)
        ax[i-1].text(140,0.45, 'n= ' +str(len(temp_units1)))
        #ax[i-1].set_title(str(Touchcount)+' TouchBlock - ' + str(Visualcount)+' VisualBloc' + str(len(temp_units1)-Touchcount-Visualcount) + ' Rest')
        #ax[i-1].grid(False)
        #ax[i-1].set_facecolor('white')
        ax[i-1].set_xticks([39.5,79.5, 119.5])
        #ax[i-1].tick_params(length=6, bottom=True)
        ax[i-1].set_xticklabels([])
        ax[i-1].vlines(39.5, 0,1)
        ax[i-1].hlines(0.5, 0,159)
        ax[i-1].set_ylim(0.4,0.7)
        #ax[i-2].set_xlabel('Seconds')
        ax[i-1].set_ylabel('mean AUC')
        ax[i-1].spines['right'].set_visible(False)
        ax[i-1].spines['top'].set_visible(False)
        #ax[i-1].legend([SomHit, VisHit, SomCR, VisCR],['SomHit','VisHit', 'SomCR', 'VisCR'])
        #ax[i-1].savefig('Mean_plot_OptoTag456_TouchStimAligned(CRs)(RespPLUS)_poster.pdf')
    fig.suptitle('DP- ALL')
    return

###############################################################################
# Figure S7C
###############################################################################

def Clusters_bargraph(mice, master_log):
    Cluster_Category=[1,2,3,4,5,6,7,8,9,10, 11, 12, 13]
    cluster_compo=[[0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0]]
    
    for i,j in enumerate(Cluster_Category):
        temp_units1=[]
        temp_units2=[]
        Touchcount=0
        Visualcount=0
        print(j)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day[0])]
                for unit in np.unique(day_log['cluster_name'][(day_log['test']==j)]):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit[0])]
                    mouse_index=[m for m,n in enumerate(mice) if n==mouse][0]
                    cluster_compo[i][mouse_index]+=1
    
    labels = ['1','2','3','4','5','6','7','8','9','10', '11', '12', '13']              
    plt.style.use('default')
    fig, ax=plt.subplots(1,1,figsize=(25,4))  
    mouse_contributions_totals=[0,0,0,0,0,0,0,0,0,0,0,0,0]    
    for i,mouse in enumerate(mice):
        mouse_contributions=[x[i] for x in cluster_compo]
        ax.bar(labels,mouse_contributions, 0.5,  bottom=mouse_contributions_totals, label=mouse)
        mouse_contributions_totals=[p+q for p,q in zip(mouse_contributions, mouse_contributions_totals)]  
        
        
    Num_per_group=[]
    Mice_in_each_group=[]    
    for Cluster_Category in [1,2,3,4,5,6,7,8,9,10, 11, 12, 13]:
        Mouse_list=[]
        Cluster_log=master_log[master_log['test']==Cluster_Category]
        for unit in np.unique(Cluster_log['unit_name']):
            unit_log=Cluster_log[Cluster_log['unit_name']==unit]
            Mouse_list.append(unit_log['mouse_name'].values[0][0])
        Num_per_group.append(np.unique(Mouse_list))
        Mice_in_each_group.append(len(np.unique(Mouse_list)))
        
        return  Num_per_group, Mice_in_each_group