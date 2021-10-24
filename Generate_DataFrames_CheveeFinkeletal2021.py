###############################################################################
# Generate_DataFrames_CheveeFinkeletal2021.py
# Proceed through this script block by block to generate:
#   - All_units_OptoMetrics_df
#   - master_log
#   - S1master_log
#   - XCorr_df
#   - S1XCorr_df
#   - master_DREADD
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



# %%
###############################################################################
# **User** update path to functions directory
###############################################################################
sys.path.append('C:\\Users\Brown Lab\Documents\Maxime_tools')

###############################################################################
################################## MATLAB #####################################
###############################################################################
########## Read Intan Data, cluster, incorporate behavior #####################
###############################################################################
#
# Run: 'EFtask_analysis_BrownLabRig_MClust.m
# Note: Ensure the function Sweep.m that is being called in the script collects 32 points per waveforms (located in +Intan>@Sweep>Line 55)
# Output: 8 .t64 files, Log, a, OptoLag, optoSpikes_log, waveform_log for each recording session
###############################################################################

# %%
###############################################################################
#Identify opto-tagged units: generates pdf plots and .txt lists
###############################################################################
from Claustrum_preprocessing_functions import Wrap_function_optoID

mouse_paths=['G:\Claustrum4', 'G:\Claustrum5', 'G:\Claustrum6']
for mouse_path in mouse_paths:
    fnall=[]
    fnall += [day for day in os.listdir(mouse_path) if day.endswith('-17')] # find all the date folders
    for day in fnall: 
    day_path=os.path.join(mouse_path, day)
    os.chdir(day_path)
    Wrap_function_optoID()
    
mouse_paths=['F:\Claustrum18', 'F:\Claustrum23',
             'F:\Claustrum25', 'G:\Claustrum31', 
             'G:\Claustrum32', 'G:\Claustrum37']
for mouse_path in mouse_paths:
    fnall=[]
    fnall += [date for date in os.listdir(mouse_path) if date.startswith('2019')] #find all the mouse folders
    for day in fnall: 
        day_path=os.path.join(mouse_path, day)
        os.chdir(day_path)
        Wrap_function_optoID()

# %%       
###############################################################################
#Identify opto-tagged units: This part is redundant with block above, here we generate a DataFrame with all metrics for all units, rather than lists of unit names
###############################################################################
# Run Generate_All_units_OptoMetrics_df.py 

# %%
###############################################################################
################################## MATLAB #####################################
###############################################################################
#Run: 'Run_Add_Trial_StartStop_to_log.m' for each mouse
###############################################################################
#Run: 'Run_OptoNetworkTagging2_SALT.m' for each mouse
###############################################################################
#RUN: 'Make_SameTTList_OptoNetwork'
###############################################################################
 
# %%
###############################################################################
#Basic initial DF construction: import data from MATLAB Log for each mouse/day
###############################################################################
columns=['mouse_name', 'date','block_type', 'trial_type', 'touch_stimulus',
 'vis_stimulus', 'response', 'trial_num', 'stim_onset', 'stim_offset',
 'licks_right', 'licks_left', 'spike_times', 'cluster_name', 'Trial_Start_Stop' ] #these are the initial columns used to build the df, more are added as more parameters are needed
master_log=pd.DataFrame(columns=columns)

mouse_path='G:\Claustrum4'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.endswith('-17')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    #log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)        


mouse_path='G:\Claustrum5'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.endswith('-17')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    #log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)        
    
mouse_path='G:\Claustrum6'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.endswith('-17')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    #log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)        

    
mouse_path='F:\Claustrum18'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
    
mouse_path='F:\Claustrum23'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
    
mouse_path='F:\Claustrum25'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
    
    
mouse_path='G:\Claustrum31'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
    
mouse_path='G:\Claustrum32'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
    
    
mouse_path='G:\Claustrum37'
fnall2=[]
fnall2 += [day for day in os.listdir(mouse_path) if day.startswith('2019')] # find all the date folders
for day in fnall2:
    day_path=os.path.join(mouse_path,day)
    day=[]
    day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
    DayLog=loadmat(day_path+'/'+str(day)[2:-2])
    log=pd.DataFrame(columns=columns)
    log['mouse_name']=DayLog['log'][:,0]
    log['date']=DayLog['log'][:,1]
    log['block_type']=DayLog['log'][:,2]
    log['trial_type']=DayLog['log'][:,3]
    log['touch_stimulus']=DayLog['log'][:,4]
    log['vis_stimulus']=DayLog['log'][:,5]
    log['response']=DayLog['log'][:,6]
    log['trial_num']=DayLog['log'][:,7]
    log['stim_onset']=DayLog['log'][:,8]
    log['stim_offset']=DayLog['log'][:,9]
    log['licks_right']=DayLog['log'][:,10]
    log['licks_left']=DayLog['log'][:,11]
    log['spike_times']=DayLog['log'][:,12]
    log['cluster_name']=DayLog['log'][:,13]
    log['Trial_Start_Stop']=DayLog['log'][:,14] #requires Add_Trial_StartStop_to_log.m to have been run
    master_log=pd.concat([master_log,log],ignore_index=True)
        
        
###############################################################################
# Add 'unit_name': indexing depends if it is Cl4/5/6 or the BrownRig mice
###############################################################################
mice=['Cl4', 'Cl5','Cl6']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for cluster in np.unique(day_log['cluster_name']):
            cluster_log=day_log.loc[np.equal(day_log['cluster_name'], cluster)]
            master_log.loc[cluster_log.index,'unit_name']= cluster_log['mouse_name'].values[0][0]+'_'+ cluster_log['date'].values[0][0]+'_'+ cluster_log['cluster_name'].values[0][0]
            print(cluster_log['mouse_name'].values[0][0]+'_'+ cluster_log['date'].values[0][0]+'_'+ cluster_log['cluster_name'].values[0][0])
        
 
mice=['Claustrum18','Claustrum23', 'Claustrum25', 'Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for cluster in np.unique(day_log['cluster_name']):
            cluster_log=day_log.loc[np.equal(day_log['cluster_name'], cluster)]
            master_log.loc[cluster_log.index,'unit_name']= cluster_log['mouse_name'].values[0][0]+'_'+ cluster_log['date'].values[0][0]+'_'+ cluster_log['cluster_name'].values[0][0]
            print(cluster_log['mouse_name'].values[0][0]+'_'+ cluster_log['date'].values[0][0]+'_'+ cluster_log['cluster_name'].values[0][0])
        

###############################################################################
#Add 'Category': OptoTag, SameTT. 
###############################################################################

#**User**: if creating the column for the first time, run this line
master_log['Category']=np.zeros(len(master_log)) #DON'T RUN WHEN JUST ADDING IN!!!! (otherwise it will erase it...)

# **User**: change the directories appropriately
#Cl4
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum4_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum4_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Cl4'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False
            
#Cl5
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum5_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum5_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Cl5'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False
#Cl6
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum6_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum6_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Cl6'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False
            
#Claustrum 18
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum18_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum18_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum18'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False

#Claustrum 23
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum23_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum23_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum23'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False     
            
#Claustrum 25
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum25_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum25_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum25'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False  
            
#Claustrum 31
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum31_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum31_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum31'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT' 
        else: 
            master_log.loc[unit_log.index,'Category']= False
            
#Claustrum 32
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum32_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum32_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum32'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False
            
#Claustrum 37
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum37_OptoList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
Units_SameTT=open('H:\\work form home 20200313\\New folder\\Claustrum37_SameTTList.txt') 
Units_SameTT= Units_SameTT.read().split('\n')
mouse='Claustrum37'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category']= 'OptoTag'
        elif unit_log['unit_name'].values[0] in Units_SameTT:
            print('SameTT')
            master_log.loc[unit_log.index,'Category']= 'SameTT'
        else: 
            master_log.loc[unit_log.index,'Category']= False
            
###############################################################################
#Add 'Category2': OptoNetwork (based on SALT test from Kepecs Lab CSHL)
###############################################################################
        
master_log['Category2']=np.zeros(len(master_log)) #DON'T RUN WHEN JUST ADDING IN!!!! (otherwise it will erase it...)

#Cl4
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum4_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl4'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        temp=unit_log['unit_name'].values[0]
        temp='Claustrum'+temp[2:] #for Cl4 Cl5 Cl6 need to do this because their unit_name is Cl4 but name in the textfile list is Claustrum4
        if temp in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
#Cl5
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum5_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl5'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        temp=unit_log['unit_name'].values[0]
        temp='Claustrum'+temp[2:] #for Cl4 Cl5 Cl6 need to do this because their unit_name is Cl4 but name in the textfile list is Claustrum4
        if temp in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False 
#Cl6
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum6_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl6'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        temp=unit_log['unit_name'].values[0]
        temp='Claustrum'+temp[2:] #for Cl4 Cl5 Cl6 need to do this because their unit_name is Cl4 but name in the textfile list is Claustrum4
        if temp in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False            
#Claustrum 18
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum18_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum18'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False


#Claustrum 23
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum23_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum23'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
        
#Claustrum 25
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum25_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum25'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
            
#Claustrum 31
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum31_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum31'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
        
#Claustrum 32
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum32_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum32'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
#Claustrum 37
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum37_OptoNetworkList_SALT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum37'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category2']= 'OptoNetwork'
        else: 
            master_log.loc[unit_log.index,'Category2']= False
            
###############################################################################
#Add 'Category3'
###############################################################################
master_log['Category3']=np.zeros(len(master_log)) #DON'T RUN WHEN JUST ADDING IN!!!! (otherwise it will erase it...)

#Cl4
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum4_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl4'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False
#Cl5
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum5_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl5'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False 
#Cl6
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum6_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl6'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False            
#Claustrum 18
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum18_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum18'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False


#Claustrum 23
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum23_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum23'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False
        
#Claustrum 25
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum25_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum25'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False
            
#Claustrum 31
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum31_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum31'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False
        
#Claustrum 32
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum32_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum32'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False
#Claustrum 37
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum37_OptoNetworkList_SALT_SameTT.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum37'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        print(unit)
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category3']= 'OptoNetwork_SameTT'
        else: 
            master_log.loc[unit_log.index,'Category3']= False            
            
###############################################################################
#Add 'Category4'
# the 'ClaustrumX_InBetweenList.txt' files were generated manually by mapping the recording
###############################################################################
master_log['Category4']=np.zeros(len(master_log)) #DON'T RUN WHEN JUST ADDING IN!!!! (otherwise it will erase it...)

#Cl4
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum4_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl4'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'     
        else: 
            master_log.loc[unit_log.index,'Category4']= False            
#Cl5
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum5_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl5'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Cl6
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum6_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Cl6'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum18
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum18_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum18'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum23
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum23_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum23'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum25
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum25_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum25'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum31
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum31_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum31'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum32
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum32_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum32'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'           
        else: 
            master_log.loc[unit_log.index,'Category4']= False
#Claustrum37
Units_OptoTag=open('H:\\work form home 20200313\\New folder\\Claustrum37_InBetweenList.txt') 
Units_OptoTag= Units_OptoTag.read().split('\n')
mouse='Claustrum37'
mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
for day in np.unique(mouse_log['date']):
    day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
    for unit in np.unique(day_log['cluster_name']):
        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
        if unit_log['unit_name'].values[0][:-5] in Units_OptoTag:
            print('OptoTag')
            master_log.loc[unit_log.index,'Category4']= 'InBetween'                       
        else: 
            master_log.loc[unit_log.index,'Category4']= False            
            
###############################################################################
# Delete all unit that are not in claustrum
###############################################################################
Master_log_OptoSameTT_LITE=master_log[(master_log['Category']=='OptoTag') | (master_log['Category']=='SameTT') | (master_log['Category2']=='OptoNetwork') | (master_log['Category3']=='OptoNetwork_SameTT')| (master_log['Category4']=='InBetween')]
master_log=Master_log_OptoSameTT_LITE        


###############################################################################
# This step is necessary for BrownLab rig mice. Because the data from the Matlab 
# logs is 1min for each trial (because I had to include the extra time to make sure
# I didn't truncate the post trial time), I had to re-apply a proper 'Stop', which
# I did in matab under the column 'Trial_Start_Stop' (Notice that was only included
# for the brownlabrig mice, not the EF mice). I create the column for the EF mice 
# also so that I can deal with just one column later on.
#
# 'licks_right_trialonly' ,'licks_left_trialonly'##
###############################################################################

#Only run the first time (initialize columns)
Name_df=['licks_right_trialonly' ,'licks_left_trialonly']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
    
#pick which mouse you want to process
mice=['Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        if i%100==0:
            print(i)
        Trial_length=master_log.loc[i,'Trial_Start_Stop'][0][1] - master_log.loc[i,'Trial_Start_Stop'][0][0] #get the length of trial Stop-start
        right_index=np.ravel(master_log.loc[i,'licks_right'][0][0]<Trial_length) #index of right licks within trial
        left_index=np.ravel(master_log.loc[i,'licks_left'][0][0]<Trial_length)#index of left licks within trial
        right_licks=master_log.loc[i,'licks_right'][0][0]
        left_licks=master_log.loc[i,'licks_left'][0][0]
        master_log.at[i,'licks_right_trialonly']=right_licks[right_index]
        master_log.at[i,'licks_left_trialonly']=left_licks[left_index]
        
mice=['Cl4','Cl5','Cl6'] #just duplicate the 'licks_right' and 'licks_left' as these don't need to be adjusted
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        if i%100==0:
            print(i)        
        master_log.at[i,'licks_right_trialonly']=  master_log.at[i,'licks_right'][0][0]
        master_log.at[i,'licks_left_trialonly']=master_log.at[i,'licks_left'][0][0]

###############################################################################
#Add 'StimALigned_spike_times'
###############################################################################       

#INIIALIZE COLUMN
Name_df=['StimALigned_spike_times']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
    
mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]  
                if not trial_log['stim_onset'][0][0].size>0:#Remove the trials where the stim was not delivered
                    master_log = master_log.drop(trial, axis=0)
                else:
                    stim_aligned_spikes=trial_log['spike_times'] - trial_log['stim_onset']
                    master_log.at[trial,'StimALigned_spike_times']=stim_aligned_spikes
                    
###############################################################################
#Add '-1to3sec_25msecbins_StimAligned'
###############################################################################       

#INIIALIZE COLUMN
Name_df=['-1to3sec_25msecbins_StimAligned']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
    
mice=['Cl4','Cl5','Cl6','Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']

for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]  
                if not trial_log['stim_onset'][0][0].size>0:#Remove the trials where the stim was not delivered 
                    master_log = master_log.drop(trial, axis=0)
                else:
                    stim_aligned_spikes=trial_log['spike_times'] - trial_log['stim_onset']
                    edges=np.arange(-1,3,0.025)
                    trialbins=np.histogram(stim_aligned_spikes, bins=edges)
                    master_log.at[trial,'-1to3sec_25msecbins_StimAligned']=trialbins[0]


################################################################################
#Build lick-aligned variables
################################################################################
         
#INIIALIZE COLUMN
Name_df=['LeftFirstLick', 'RightFirstLick', 'FirstLick']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing

mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]
                       
                if trial_log['stim_onset']: 
            
                    if not trial_log['licks_left_trialonly'].size: #if there was no lick, make it 0 such that it'll be sorted first
                        master_log.at[trial,'LeftFirstLick']=0; # I can use trial as the index directly to master_log 
                    else:
                        TempLicks=trial_log['licks_left_trialonly']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                        #TempLicks=TempLicks[0][0] #get rid of extra brackets
                        if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                            master_log.at[trial,'LeftFirstLick']=0;
                        else:
                            master_log.at[trial,'LeftFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count

                    if not trial_log['licks_right_trialonly'].size: #if there was no lick, make it 0 such that it'll be sorted first
                        master_log.at[trial,'RightFirstLick']=0; # I can use trial as the index directly to master_log 
                    else:
                        TempLicks=trial_log['licks_right_trialonly']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                        #TempLicks=TempLicks[0][0] #get rid of extra brackets
                        if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                            master_log.at[trial,'RightFirstLick']=0;
                        else:
                            master_log.at[trial,'RightFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count
#The following loop goes through the possibilites and populates 
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]
                
                if trial_log['stim_onset']: 
        
                    if (trial_log['LeftFirstLick']!=0) | (trial_log['RightFirstLick']!=0):
                        if (trial_log['LeftFirstLick']==0) & (trial_log['RightFirstLick']!=0):
                            master_log.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                        elif (trial_log['LeftFirstLick']!=0) & (trial_log['RightFirstLick']==0):
                            master_log.at[trial, 'FirstLick']=trial_log['LeftFirstLick']
                        elif trial_log['LeftFirstLick']> trial_log['RightFirstLick']:
                            master_log.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                        elif trial_log['LeftFirstLick']< trial_log['RightFirstLick']:
                            master_log.at[trial, 'FirstLick']=trial_log['LeftFirstLick']
                            
###############################################################################
#Add '-2to2sec_25msecbins_LickAligned'
############################################################################### 
        
#INIIALIZE COLUMN
Name_df=['-2to2sec_25msecbins_LickAligned', 'LickALigned_spike_times']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end preprocessing

mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial] 
                lick_aligned_spikes=master_log.loc[trial,'StimALigned_spike_times'] - master_log.loc[trial,'FirstLick']
                edges=np.arange(-2,2,0.025)
                trialbins,_=np.histogram(lick_aligned_spikes, bins=edges)
                master_log.at[trial,'LickALigned_spike_times']=lick_aligned_spikes
                master_log.at[trial,'-2to2sec_25msecbins_LickAligned']=trialbins

###############################################################################
#Add ''Zscored_-1to3sec_25msecbins_StimAligned''
############################################################################### 
        
#INIIALIZE COLUMN
Name_df=['Zscored_-1to3sec_25msecbins_StimAligned']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing       

mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit) 
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            
            #get the '-1to3sec_25msecbins_StimAligned' from each trial for the unit and Zscore
            temp=np.array(0); #initialize an array (this 0 will be extra but inconsequential)
            for each in unit_log['-1to3sec_25msecbins_StimAligned']:
                unit_index=unit_log.index
                temp=np.hstack((temp, each)) 
            #plt.scatter(np.arange(len(temp)), temp)
            Zscored_temp=sp.stats.zscore(temp)
            #Now repackage the array into a list in the original order
            Zscored_bins=[]; #initialize list
            for m in np.arange(1,len(unit_log)*159,159):
                Zscored_bins.append(Zscored_temp[m:m+159])
            #put it back into master_log at the right index
            for i,index in enumerate(unit_index):
                master_log.at[index, 'Zscored_-1to3sec_25msecbins_StimAligned']=Zscored_bins[i]
                
                
###############################################################################
#Add ''Zscored_-2to2sec_25msecbins_LickAligned'
############################################################################### 
        
#INIIALIZE COLUMN
Name_df=['Zscored_-2to2sec_25msecbins_LickAligned']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing       

mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit) 
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            
            #get the '-1to3sec_25msecbins_StimAligned' from each trial for the unit and Zscore
            temp=np.array(0); #initialize an array (this 0 will be extra but inconsequential)
            for each in unit_log['-2to2sec_25msecbins_LickAligned']:
                unit_index=unit_log.index
                temp=np.hstack((temp, each)) 
            #plt.scatter(np.arange(len(temp)), temp)
            Zscored_temp=sp.stats.zscore(temp)
            #Now repackage the array into a list in the original order
            Zscored_bins=[]; #initialize list
            for m in np.arange(1,len(unit_log)*159,159):
                Zscored_bins.append(Zscored_temp[m:m+159])
            #put it back into master_log at the right index
            for i,index in enumerate(unit_index):
                master_log.at[index, 'Zscored_-2to2sec_25msecbins_LickAligned']=Zscored_bins[i]
                
                               
###############################################################################
#Add 'Stim/Block/Response'
############################################################################### 
#For the Original task: Som:lick right, vis:lick left       
#For the REVERSED task: Som:lick left, vis:lick right
mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37'] 
for mouse in mice:  
    for index in master_log[(master_log['mouse_name']==mouse)].index:
        if index%100==0:
            print(index) 
        if [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index, 'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==1)][0]:        
            master_log.loc[index, 'Stim/Block/Response']='SomHit'
        elif [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomMiss'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==2)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisHit'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisMiss'
        elif [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomCR'
        elif [(( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==2))
                    | (( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==1))
                    | (( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==2))][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomFA'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisCR'
        elif [(( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==2))
                    | (( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==1))
                    | (( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==1))][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisFA'

mice=['Claustrum18', 'Claustrum23', 'Claustrum25'] 
for mouse in mice:  
    for index in master_log[(master_log['mouse_name']==mouse)].index:
        if index%100==0:
            print(index)
        if [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index, 'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==2)][0]:        
            master_log.loc[index, 'Stim/Block/Response']='SomHit'
        elif [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomMiss'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==1)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisHit'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisMiss'
        elif [(master_log.loc[index,'block_type']=='Whisker') & (master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomCR'
        elif [(( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==2))
                    | (( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==1))
                    | (( master_log.loc[index,'block_type']=='Whisker') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==1))][0]:
            master_log.loc[index, 'Stim/Block/Response']='SomFA'
        elif [(master_log.loc[index,'block_type']=='Visual') & (master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_log.loc[index,'response']==0)][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisCR'
        elif [(( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==2))
                    | (( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_log.loc[index,'response']==1))
                    | (( master_log.loc[index,'block_type']=='Visual') & ( master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_log.loc[index,'response']==2))][0]:
            master_log.loc[index, 'Stim/Block/Response']='VisFA'
            
               
###############################################################################
#'Reward'
###############################################################################       
mice=['Cl4','Cl5','Cl6', 'Claustrum18','Claustrum23','Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:                
    for index in master_log[(master_log['mouse_name']==mouse)].index:
        if index%100==0:
            print(index)
        if [(master_log.loc[index,'Stim/Block/Response']=='SomHit')][0] | [(master_log.loc[index, 'Stim/Block/Response']=='VisHit')][0]:
             master_log.loc[index, 'Reward']=1
        else:
            master_log.loc[index, 'Reward']=0
            
            


###############################################################################
# 'rewarded_licks','non_rewarded_licks'
###############################################################################           
        
#'rewarded_licks' and 'unrewarded_licks'
Name_df=['rewarded_licks','non_rewarded_licks']
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
                

mice=['Cl4','Cl5','Cl6'] 
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        if master_log.loc[i,'Reward']==1:
            if master_log.loc[i,'block_type']=='Visual':
                a=(master_log.loc[i,'licks_left_trialonly'] - master_log.loc[i,'stim_onset'])>0.01
                b=(master_log.loc[i,'licks_left_trialonly'] - master_log.loc[i,'stim_onset'])<3
                mask = [all(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                temp1=np.ravel(master_log.loc[i,'licks_left_trialonly'][(master_log.loc[i,'licks_left_trialonly'] - master_log.loc[i,'stim_onset'])<0.01])# licks from that trial before the stim are unrewarded
                temp2=np.ravel(master_log.loc[i,'licks_right_trialonly']) # licks on wrong port
                non_rewarded_licks=np.hstack((temp1, temp2))
                master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
            elif master_log.loc[i,'block_type']=='Whisker':
                a=(master_log.loc[i,'licks_right_trialonly'] - master_log.loc[i,'stim_onset'])>0.01
                b=(master_log.loc[i,'licks_right_trialonly'] - master_log.loc[i,'stim_onset'])<3
                mask = [all(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'][(master_log.loc[i,'licks_right_trialonly'] - master_log.loc[i,'stim_onset'])<0.01])# licks from that trial before the stim are unrewarded
                temp2=np.ravel(master_log.loc[i,'licks_left_trialonly']) # licks on wrong port
                non_rewarded_licks=np.hstack((temp1, temp2))
                master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
        else:
            master_log.at[i,'rewarded_licks']=[]
            temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'])
            temp2=np.ravel(master_log.loc[i,'licks_left_trialonly'])
            non_rewarded_licks=np.hstack((temp1, temp2))
            master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...

mice=['Claustrum31', 'Claustrum32', 'Claustrum37'] 
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit) 
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for j,i in enumerate(unit_log.index):
            
                if master_log.loc[i,'Reward']==1:
                    if master_log.loc[i,'block_type']=='Visual':
                        a=(master_log.loc[i,'licks_left_trialonly'] - master_log.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1:
                            b=(master_log.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_log.loc[i,'block_type']=='Whisker':
                        a=(master_log.loc[i,'licks_right_trialonly'] - master_log.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: 
                            b=(master_log.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                else: #no reward
                    if master_log.loc[i,'block_type']=='Visual':
                        a=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        if j==0:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: 
                            b=(master_log.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_log.loc[i,'block_type']=='Whisker':
                        a=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        if j==0:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: 
                            b=(master_log.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
              
                            
mice=[ 'Claustrum18', 'Claustrum23', 'Claustrum25'] 
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit) 
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for j,i in enumerate(unit_log.index):
            
                if master_log.loc[i,'Reward']==1:
                    if master_log.loc[i,'block_type']=='Whisker':
                        a=(master_log.loc[i,'licks_left_trialonly'] - master_log.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: #This stratgey very slightly overcounts rewarded licks, but it's ok given the analysis we will use this for
                            b=(master_log.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_log.loc[i,'block_type']=='Visual':
                        a=(master_log.loc[i,'licks_right_trialonly'] - master_log.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: #This stratgey very slightly overcounts rewarded licks, but it's ok given the analysis we will use this for
                            b=(master_log.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                else: #no reward
                    if master_log.loc[i,'block_type']=='Whisker':
                        a=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        if j==0:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: #This stratgey very slightly overcounts rewarded licks, but it's ok given the analysis we will use this for
                            b=(master_log.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_log.loc[i,'block_type']=='Visual':
                        a=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        if j==0:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_log.loc[unit_log.index[j-1],'Reward']==1: #This stratgey very slightly overcounts rewarded licks, but it's ok given the analysis we will use this for
                            b=(master_log.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_log.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_log.at[i,'rewarded_licks']=master_log.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_log.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_log.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
              
###############################################################################
#'non_rewarded_Right_licks', 'non_rewarded_Left_licks'
############################################################################### 

#'non_rewarded_Right_licks', 'non_rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licks', 'non_rewarded_Left_licks'    ]
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
    
mice=['Cl4','Cl5','Cl6', 'Claustrum18', 'Claustrum23', 'Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37']
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        temp_right=[]
        temp_left=[]
        trial_log=master_log.loc[i,:]
        for ind in range(0, len(trial_log.non_rewarded_licks)):
            
            if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                temp_right.append(trial_log['non_rewarded_licks'][ind])
       
            if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                temp_left.append(trial_log['non_rewarded_licks'][ind])
        master_log.at[i,'non_rewarded_Right_licks']=temp_right
        master_log.at[i,'non_rewarded_Left_licks']=temp_left

###############################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'  
############################################################################### 
            
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['rewarded_Right_licks', 'rewarded_Left_licks'    ]
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
#end of preprocessing
mice=['Cl4','Cl5','Cl6','Claustrum18', 'Claustrum23', 'Claustrum25','Claustrum31', 'Claustrum32', 'Claustrum37'] 
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        temp_right=[]
        temp_left=[]
        trial_log=master_log.loc[i,:]
        for ind in range(0, len(trial_log.rewarded_licks)):
            
            if trial_log['rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                temp_right.append(trial_log['rewarded_licks'][ind])
       
            if trial_log['rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                temp_left.append(trial_log['rewarded_licks'][ind])
        master_log.at[i,'rewarded_Right_licks']=temp_right
        master_log.at[i,'rewarded_Left_licks']=temp_left
        
###############################################################################
        # Strategy for creating columns in blocks below is different #
###############################################################################
        # 'NR_licks_aligned_spikes', 'Rewarded_licks_aligned_spikes'
        # 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
        # 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
############################################################################### 
###############################################################################
        
# 'NR_licks_aligned_spikes'
#Note: only the licks with 'First' in their names are relative to stim onset. All others (as in block below) are relative to trial start
new = []
for num in master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    new_col = []
    for ind in range(0, len(master_log.non_rewarded_licks[num])):
        new_col.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))
    new.append(new_col)
master_log['NR_licks_aligned_spikes'] = new

# 'Rewarded_licks_aligned_spikes'
new = []
for num in master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    new_col = []
    for ind in range(0, len(master_log.rewarded_licks[num])):
        new_col.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['rewarded_licks'][num][ind]]))
    new.append(new_col)
master_log['Rewarded_licks_aligned_spikes'] = new

# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'     
NR_Right_licks_aligned_spikes = []
NR_Left_licks_aligned_spikes = []

for num in master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    temp_right = []
    temp_left = []
    for ind in range(0, len(master_log.non_rewarded_licks[num])):
        if master_log['non_rewarded_licks'][num][ind] in  master_log['licks_right'][num][0][0]:
            temp_right.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))
        elif master_log['non_rewarded_licks'][num][ind] in  master_log['licks_left'][num][0][0]:
            temp_left.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))
        else:
            print('Abandoned?')
                
    NR_Right_licks_aligned_spikes.append(temp_right)
    NR_Left_licks_aligned_spikes.append(temp_left)
master_log['NR_Right_licks_aligned_spikes'] = NR_Right_licks_aligned_spikes
master_log['NR_Left_licks_aligned_spikes'] = NR_Left_licks_aligned_spikes

# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
Rewarded_Right_licks_aligned_spikes = []
Rewarded_Left_licks_aligned_spikes = []

for num in master_log.rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    temp_right = []
    temp_left = []
    for ind in range(0, len(master_log.rewarded_licks[num])):
        if master_log['rewarded_licks'][num][ind] in  master_log['licks_right'][num][0][0]:
            temp_right.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['rewarded_licks'][num][ind]]))
        elif master_log['rewarded_licks'][num][ind] in  master_log['licks_left'][num][0][0]:
            temp_left.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['rewarded_licks'][num][ind]]))
        else:
            print('Abandoned?')
                
    Rewarded_Right_licks_aligned_spikes.append(temp_right)
    Rewarded_Left_licks_aligned_spikes.append(temp_left)
master_log['Rewarded_Right_licks_aligned_spikes'] = Rewarded_Right_licks_aligned_spikes
master_log['Rewarded_Left_licks_aligned_spikes'] = Rewarded_Left_licks_aligned_spikes


##########################################################################################################
## Repeat the non_rewarded lick analysis but only count the licks that are not in the middle of bouts. So only isolated licks or first lick of a bout
##########################################################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licksNonBouts', 'non_rewarded_Left_licksNonBouts'    ]
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
    
mice=[ 'Claustrum18','Claustrum23', 'Claustrum25', 'Cl4','Cl5', 'Cl6','Claustrum31','Claustrum32', 'Claustrum37']
for mouse in mice:
    for i in master_log[master_log['mouse_name']==mouse].index:
        temp_right=[]
        temp_left=[]
        trial_log=master_log.loc[i,:]
        for ind in range(0, len(trial_log.non_rewarded_licks)):
            if ind==0: #if first lick analyzed no need to check before
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                    temp_right.append(trial_log['non_rewarded_licks'][ind])
               
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                    temp_left.append(trial_log['non_rewarded_licks'][ind])
            elif trial_log['non_rewarded_licks'][ind] - trial_log['non_rewarded_licks'][ind-1]>0.5: # so only if the lick is not preceded by a lick within 500msec
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                    temp_right.append(trial_log['non_rewarded_licks'][ind])
               
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                    temp_left.append(trial_log['non_rewarded_licks'][ind])
        master_log.at[i,'non_rewarded_Right_licksNonBouts']=temp_right
        master_log.at[i,'non_rewarded_Left_licksNonBouts']=temp_left
        test=master_log['non_rewarded_Left_licksNonBouts']
        
# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
NR_Right_licks_aligned_spikes_NonBouts = []
NR_Left_licks_aligned_spikes_NonBouts = []

for num in master_log.non_rewarded_licks.index:
    temp_right = []
    temp_left = []
    for ind in range(0, len(master_log.non_rewarded_licks[num])):
        if master_log['non_rewarded_licks'][num][ind] in  master_log['non_rewarded_Right_licksNonBouts'][num]:
            temp_right.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))
        elif master_log['non_rewarded_licks'][num][ind] in  master_log['non_rewarded_Left_licksNonBouts'][num]:
            temp_left.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))

                
    NR_Right_licks_aligned_spikes_NonBouts.append(temp_right)
    NR_Left_licks_aligned_spikes_NonBouts.append(temp_left)
master_log['NR_Right_licks_aligned_spikes_NonBouts'] = NR_Right_licks_aligned_spikes_NonBouts
master_log['NR_Left_licks_aligned_spikes_NonBouts'] = NR_Left_licks_aligned_spikes_NonBouts

##########################################################################################################
## Repeat the ITBNs but adding the limitation of no lick "Following" on top of no lick "preceding"
##########################################################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licksNonBouts2', 'non_rewarded_Left_licksNonBouts2'    ]
for i,name in enumerate(Name_df):   
    master_log[Name_df[i]]=np.zeros(len(master_log))
    master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
    
mice=[ 'Claustrum18','Claustrum23', 'Claustrum25', 'Cl4','Cl5', 'Cl6','Claustrum31','Claustrum32', 'Claustrum37']
for mouse in mice:
    print(mouse)
    for i in master_log[master_log['mouse_name']==mouse].index:
        temp_right=[]
        temp_left=[]
        trial_log=master_log.loc[i,:]
        if len(trial_log.non_rewarded_licks)==1: #if there are only 1 lick, count it
            if trial_log['non_rewarded_licks'][0] in  trial_log['licks_right'][0][0]:
                    temp_right.append(trial_log['non_rewarded_licks'][0])
            if trial_log['non_rewarded_licks'][0] in  trial_log['licks_left'][0][0]:
                    temp_left.append(trial_log['non_rewarded_licks'][0])
        if len(trial_log.non_rewarded_licks)==2: #if there are only 2 lick, count them only if sep by 500msec
            if (trial_log['non_rewarded_licks'][1] -  trial_log['non_rewarded_licks'][0])>0.5:
                for ind in range(0, len(trial_log.non_rewarded_licks)):
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                            temp_right.append(trial_log['non_rewarded_licks'][ind])
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                            temp_left.append(trial_log['non_rewarded_licks'][ind])
        else:
            for ind in range(0, len(trial_log.non_rewarded_licks))[1:-1]:
                if ((trial_log['non_rewarded_licks'][ind] - trial_log['non_rewarded_licks'][ind-1]>0.5) & (trial_log['non_rewarded_licks'][ind+1] - trial_log['non_rewarded_licks'][ind]>0.5)): # so only if the lick isolated by 500msec
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                        temp_right.append(trial_log['non_rewarded_licks'][ind])
                       
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                        temp_left.append(trial_log['non_rewarded_licks'][ind])
                    
        master_log.at[i,'non_rewarded_Right_licksNonBouts2']=temp_right
        master_log.at[i,'non_rewarded_Left_licksNonBouts2']=temp_left
        
# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
NR_Right_licks_aligned_spikes_NonBouts = []
NR_Left_licks_aligned_spikes_NonBouts = []

for num in master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    temp_right = []
    temp_left = []
    for ind in range(0, len(master_log.non_rewarded_licks[num])):
        if master_log['non_rewarded_licks'][num][ind] in  master_log['non_rewarded_Right_licksNonBouts2'][num]:
            temp_right.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))
        elif master_log['non_rewarded_licks'][num][ind] in  master_log['non_rewarded_Left_licksNonBouts2'][num]:
            temp_left.append(np.array([np.squeeze(master_log['spike_times'][num])-master_log['non_rewarded_licks'][num][ind]]))

                
    NR_Right_licks_aligned_spikes_NonBouts.append(temp_right)
    NR_Left_licks_aligned_spikes_NonBouts.append(temp_left)
master_log['NR_Right_licks_aligned_spikes_NonBouts2'] = NR_Right_licks_aligned_spikes_NonBouts
master_log['NR_Left_licks_aligned_spikes_NonBouts2'] = NR_Left_licks_aligned_spikes_NonBouts

##########################################################################################################
## Add 'LickPref'
##########################################################################################################

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
            


master_counter=0
master_LickPref=[]
Unit_names=[]
#trial_types=['SomHit','VisHit','SomFA','VisFA','SomCR','VisCR','SomMiss','VisMiss']     [ 'VisHit', 'VisFA','SomHit', 'SomFA']   ['SomMiss','VisMiss','SomCR','VisCR'] 

for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']): 
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            Unit_names.append(unit_log['unit_name'].values[0])
            
            print(unit_log['unit_name'].values[0])

            master_LickPref.append((peak_Right[master_counter]-peak_Left[master_counter])/(peak_Right[master_counter]+peak_Left[master_counter]))
            master_counter+=1

#Add lickPref to master_log
count=0
for mouse in mice:
    mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            master_log.at[unit_log.index, 'LickPref']=master_LickPref[count]
            count+=1      


###############################################################################
# Clean up
###############################################################################
########################
#DROP
########################

#ISI violations
units_to_drop=['Claustrum18_20190402_TT4clst3',
 'Claustrum37_20191123_TT2clst3', #OptoTag
 'Cl4_05-23-17_TT5clst4',
 'Cl4_05-26-17_TT2clst4', #OptoTag
 'Cl4_05-26-17_TT6clst3', #SALT
 'Cl4_06-03-17_TT4clst2',
 'Cl4_06-08-17_TT7clst2',
 'Cl4_06-09-17_TT5clst1', #OptoTag
 'Cl5_06-09-17_TT6clst1',
 'Claustrum18_20190401_TT6clst1', #OptoTag
 'Claustrum18_20190420_TT7clst3',
 'Claustrum18_20190421_TT5clst2', #SALT
 'Claustrum23_20190726_TT7clst1', #OptoTag
 'Claustrum23_20190730_TT3clst1',
 'Claustrum25_20190801_TT7clst2',
 'Claustrum25_20190804_TT7clst2', #SALT
 'Claustrum37_20191117_TT3clst3', #OptoTag
 'Claustrum37_20191122_TT3clst2', #OptoTag
 'Claustrum37_20191125_TT2clst2',
 'Cl4_05-26-17_TT2clst1',
 'Cl4_05-26-17_TT2clst2',
 'Cl4_05-26-17_TT2clst3', #SALT
 'Cl4_06-08-17_TT7clst1',
 'Cl4_06-08-17_TT7clst3',
 'Claustrum18_20190402_TT6clst1',
 'Claustrum18_20190402_TT6clst2',
 'Claustrum18_20190402_TT6clst3',
 'Claustrum37_20191117_TT3clst1',
 'Claustrum37_20191117_TT3clst2',
 'Claustrum37_20191119_TT3clst1',
 'Claustrum37_20191120_TT3clst1',
 'Claustrum37_20191120_TT3clst2',
 'Claustrum37_20191120_TT3clst3']
  
index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)

########################
#Funky waveform or drift
########################
units_to_drop=['Cl4_05-24-17_TT8clst2',
'Cl4_06-02-17_TT4clst2', #OptoTag
'Claustrum25_20190801_TT2clst1', #SALT
'Claustrum31_20191114_TT8clst2',
'Claustrum31_20191105_TT3clst2',
'Claustrum31_20191105_TT3clst3',
'Claustrum31_20191106_TT4clst3',
'Claustrum31_20191106_TT4clst4',
'Claustrum37_20191121_TT8clst4',
'Cl4_05-19-17_TT2clst1',
'Cl4_05-24-17_TT1clst3',
'Cl4_05-24-17_TT2clst3',
'Cl4_05-25-17_TT2clst3',
'Cl4_06-08-17_TT1clst1',
'Cl4_06-08-17_TT1clst4',
'Cl4_06-08-17_TT1clst6',
'Cl5_05-25-17_TT2clst4',
'Cl5_05-25-17_TT2clst5',
'Cl5_05-25-17_TT2clst6',
'Cl5_05-25-17_TT5clst1', #SALT
'Cl5_05-25-17_TT8clst2',
'Claustrum18_20190423_TT1clst5',
'Claustrum23_20190725_TT8clst2',
'Claustrum23_20190726_TT8clst2',
'Claustrum23_20190726_TT4clst3',]

index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)

########################
#Lick artifact
########################
units_to_drop=['Claustrum25_20190806_TT1clst4', #OptoTag
               'Claustrum25_20190812_TT1clst3'] #OptoTag

index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)

########################
#Duplicate by CCG
########################
units_to_drop=['Cl4_05-19-17_TT5clst2', #clustering error: same spikes belong to clst2 and clst3
               'Cl4_05-19-17_TT5clst3', #clustering error: same spikes belong to clst2 and clst3
               'Cl4_05-26-17_TT6clst5', #clustering error: same unit as TT6clst4
               'Claustrum23_20190730_TT6clst2', #clustering error: same as clst1 (there was a typo:20190703 instead of 20190730, so it is not dropped in dataset, drop it)
               'Claustrum31_20191108_TT5clst3', #Not in claustrum
               'Cl4_05-23-17_TT1clst1', #clustering error: same as clst3
               'Cl5_06-04-17_TT2clst4',  #clustering error: same as clst5
               'Cl4_05-23-17_TT1clst1',#clustering error: same as clst3
               'Cl5_06-04-17_TT2clst4'] #clustering error: same as clst5
index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)

######################
# Removed consequences of other removed
######################
units_to_drop=[
'Cl4_06-09-17_TT5clst3',
'Cl4_05-26-17_TT2clst2',
 'Claustrum25_20190801_TT2clst2',
 'Claustrum37_20191117_TT3clst1',
 'Claustrum37_20191117_TT3clst2']
index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)

units_to_drop=[]
for unit in np.unique(master_log['unit_name']):
    if unit[:-5]=='Claustrum18_20190402_TT6':
        units_to_drop.append(unit)
    elif unit[:-5]=='Claustrum25_20190804_TT2':
        units_to_drop.append(unit)
index_to_drop=[]
for each in units_to_drop:
  unit_log=master_log[master_log['unit_name']==each]
  index_to_drop.append(unit_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
master_log.drop(indices, axis=0, inplace=True)


##########################################################################################################
## Add AUC results to master_log
##########################################################################################################
from Claustrum_preprocessing_functions import Add_AUC_to_master_log
Add_AUC_to_master_log(master_log)


##########################################################################################################
## Cluster response and add group ID under column name 'test'
##########################################################################################################
from Claustrum_preprocessing_functions import clustering_add_test
clustering_add_test(master_log)


##########################################################################################################
## SAVE
##########################################################################################################

master1a=master_log.loc[master_log.index[0:20000],:]
master1a.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partA.pkl', protocol=4)
del master1a

master1b=master_log.loc[master_log.index[20000:40000],:]
master1b.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partB.pkl', protocol=4)
del master1b

master1c=master_log.loc[master_log.index[40000:60000],:]
master1c.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partC.pkl', protocol=4)
del master1c

master1d=master_log.loc[master_log.index[60000:80000],:]
master1d.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partD.pkl', protocol=4)
del master1d

master1e=master_log.loc[master_log.index[80000:100000],:]
master1e.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partE.pkl', protocol=4)
del master1e

master1f=master_log.loc[master_log.index[100000:120000],:]
master1f.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partF.pkl', protocol=4)
del master1f

master1g=master_log.loc[master_log.index[120000:140000],:]
master1g.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partG.pkl', protocol=4)
del master1g

master1h=master_log.loc[master_log.index[140000:160000],:]
master1h.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partH.pkl', protocol=4)
del master1h

master1i=master_log.loc[master_log.index[160000:180000],:]
master1i.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partI.pkl', protocol=4)
del master1i

master1j=master_log.loc[master_log.index[180000:200000],:]
master1j.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partJ.pkl', protocol=4)
del master1j

master1k=master_log.loc[master_log.index[200000:],:]
master1k.to_pickle( 'C:\\Users\\Brown Lab\\Desktop\\Maxime_ClaustrumRevision\\test\\master_log-partK.pkl', protocol=4)
del master1k

# %%
##############################################################################
##############################################################################
# Process S1 data from Eric Finkel
##############################################################################
##############################################################################

# Data frame provided by EF: S1log_df.h5
S1master_log = pd.read_hdf('C:\\Users\Brown Lab\Desktop\Analysis working folder\Claustrum explore\S1log_df.h5','table')
mice_to_drop=[ 'EF0081', 'EF0083',
       'EF0084', 'EF0085', 'EF0088', 'EF0089', 'EF0091', 'EF0099',
       'EF0101', 'EF0102', 'EF0111', 'EF0112', 'EF0114'] #Mice from EF paper doing other things
index_to_drop=[]
for each in mice_to_drop:
  mouse_log=S1master_log[S1master_log['mouse_name']==each]
  index_to_drop.append(mouse_log.index.values)
  
indices=[x for i in index_to_drop for x in i]
S1master_log.drop(indices, axis=0, inplace=True)


#Create the columns in master_log to be filled  (ONLY RUN ONCE! - which you've done already... rerunning will erase what currently exists)
Name_df=['-1to3sec_25msecbins_StimAligned']
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
#end of preprocessing

List_to_Del=[] #we save the index of trials with no stim and delete them at the end
#ADD DATA
for mouse in np.unique(S1master_log['mouse_name']):
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]  
                if trial_log['stim_onset']==0:#Remove the trials where the stim was not delivered
                    List_to_Del.append(trial)
                else:
                    stim_aligned_spikes=trial_log['spike_times'] - trial_log['stim_onset']
                    edges=np.arange(-1,3,0.025)
                    trialbins=np.histogram(stim_aligned_spikes, bins=edges)
                    S1master_log.at[trial,'-1to3sec_25msecbins_StimAligned']=trialbins[0]

S1master_log = S1master_log.drop(List_to_Del, axis=0)

######################################################################################################################
#'Stim/Block/Response'
######################################################################################################################
#### original task: Som:lick rigt, Vis:lick left

mice=['EF0074', 'EF0076','EF0077','EF0079'] 
for mouse in mice:
    for index in S1master_log[S1master_log['mouse_name']==mouse].index:
        if index%100==0:
            print(index)
        if [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index, 'trial_type']=='Stim_Som_NoCue') & (S1master_log.loc[index,'response']==1)][0]:        
            S1master_log.loc[index, 'Stim/Block/Response']='SomHit'
        elif [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomMiss'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (S1master_log.loc[index,'response']==2)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisHit'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisMiss'
        elif [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomCR'
        elif [(( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( S1master_log.loc[index,'response']==1))
                    | (( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( S1master_log.loc[index,'response']==2))
                    | (( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( S1master_log.loc[index,'response']==2))][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomFA'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisCR'
        elif [(( S1master_log.loc[index,'block_type']=='Visal') & ( S1master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( S1master_log.loc[index,'response']==1))
                    | (( S1master_log.loc[index,'block_type']=='Visual') & ( S1master_log.loc[index,'trial_type']=='Stim_Som_NoCue') & ( S1master_log.loc[index,'response']==2))
                    | (( S1master_log.loc[index,'block_type']=='Visual') & ( S1master_log.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( S1master_log.loc[index,'response']==1))][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisFA'
            
#Some of the trial_types are  '1CycStim_Som_NoCue', '1CycStim_Vis_NoCue', 'Stim_Som_NoCue', 'Stim_Vis_NoCue'
mice=['EF0074', 'EF0076','EF0077','EF0079']
for mouse in mice:
    for index in S1master_log[S1master_log['mouse_name']==mouse].index:
        if index%100==0:
            print(index)
        if [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index, 'trial_type']=='1CycStim_Som_NoCue') & (S1master_log.loc[index,'response']==1)][0]:        
            S1master_log.loc[index, 'Stim/Block/Response']='SomHit'
        elif [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index,'trial_type']=='1CycStim_Som_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomMiss'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & (S1master_log.loc[index,'response']==2)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisHit'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisMiss'
        elif [(S1master_log.loc[index,'block_type']=='Whisker') & (S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomCR'
        elif [(( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & ( S1master_log.loc[index,'response']==1))
                    | (( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & ( S1master_log.loc[index,'response']==2))
                    | (( S1master_log.loc[index,'block_type']=='Whisker') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Som_NoCue') & ( S1master_log.loc[index,'response']==2))][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='SomFA'
        elif [(S1master_log.loc[index,'block_type']=='Visual') & (S1master_log.loc[index,'trial_type']=='1CycStim_Som_NoCue') & (S1master_log.loc[index,'response']==0)][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisCR'
        elif [(( S1master_log.loc[index,'block_type']=='Visal') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Som_NoCue') & ( S1master_log.loc[index,'response']==1))
                    | (( S1master_log.loc[index,'block_type']=='Visual') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Som_NoCue') & ( S1master_log.loc[index,'response']==2))
                    | (( S1master_log.loc[index,'block_type']=='Visual') & ( S1master_log.loc[index,'trial_type']=='1CycStim_Vis_NoCue') & ( S1master_log.loc[index,'response']==1))][0]:
            S1master_log.loc[index, 'Stim/Block/Response']='VisFA'
 

#################################################################################
## Add 'FirstLick' Column
#################################################################################
Name_df=['LeftFirstLick',
                     'RightFirstLick',
                     'FirstLick']
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
#end of preprocessing
for mouse in mice:
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]
                       
                if trial_log['stim_onset']: 
            
                    if not trial_log['licks_left'].size: #if there was no lick, make it 0 such that it'll be sorted first
                        S1master_log.at[trial,'LeftFirstLick']=0; # I can use trial as the index directly to master_log 
                    else:
                        TempLicks=trial_log['licks_left']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                        if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                            S1master_log.at[trial,'LeftFirstLick']=0;
                        else:
                            S1master_log.at[trial,'LeftFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count

                    if not trial_log['licks_right'].size: #if there was no lick, make it 0 such that it'll be sorted first
                        S1master_log.at[trial,'RightFirstLick']=0; # I can use trial as the index directly to master_log 
                    else:
                        TempLicks=trial_log['licks_right']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                        if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                            S1master_log.at[trial,'RightFirstLick']=0;
                        else:
                            S1master_log.at[trial,'RightFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count


#The followinf loop goes through the possibilites and populates the
for mouse in mice:
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]
                
                if trial_log['stim_onset']:
        
                    if (trial_log['LeftFirstLick']!=0) | (trial_log['RightFirstLick']!=0):
                        if (trial_log['LeftFirstLick']==0) & (trial_log['RightFirstLick']!=0):
                            S1master_log.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                        elif (trial_log['LeftFirstLick']!=0) & (trial_log['RightFirstLick']==0):
                            S1master_log.at[trial, 'FirstLick']=trial_log['LeftFirstLick']
                        elif trial_log['LeftFirstLick']> trial_log['RightFirstLick']:
                            S1master_log.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                        elif trial_log['LeftFirstLick']< trial_log['RightFirstLick']:
                            S1master_log.at[trial, 'FirstLick']=trial_log['LeftFirstLick']

###############################################################################
        #Add 'StimALigned_spike_times'
############################################################################### 
#INIIALIZE COLUMN
Name_df=['StimALigned_spike_times']
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
#end of preprocessing
    
for mouse in mice:
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial]  
                if not trial_log['stim_onset'].size>0:
                    S1master_log = S1master_log.drop(trial, axis=0)
                else:
                    stim_aligned_spikes=trial_log['spike_times'] - trial_log['stim_onset']
                    S1master_log.at[trial,'StimALigned_spike_times']=stim_aligned_spikes
                    
###############################################################################
#Add '-2to2sec_25msecbins_LickAligned' and 'LickALigned_spike_times'
############################################################################### 
        
#INIIALIZE COLUMN
Name_df=['-2to2sec_25msecbins_LickAligned', 'LickALigned_spike_times']
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
#end preprocessing

for mouse in mice:
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for unit in np.unique(day_log['cluster_name']):
            print(unit)
            unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
            for trial in unit_log.index:
                trial_log=unit_log.loc[trial] 
                lick_aligned_spikes=S1master_log.loc[trial,'StimALigned_spike_times'] - S1master_log.loc[trial,'FirstLick']
                edges=np.arange(-2,2,0.025)
                trialbins,_=np.histogram(lick_aligned_spikes, bins=edges)
                S1master_log.at[trial,'LickALigned_spike_times']=lick_aligned_spikes
                S1master_log.at[trial,'-2to2sec_25msecbins_LickAligned']=trialbins
                

###############################################################################
# 'Reawrd'
###############################################################################    
for mouse in mice:                
    for index in S1master_log[(S1master_log['mouse_name']==mouse)].index:
        if index%100==0:
            print(index)
        if [(S1master_log.loc[index,'Stim/Block/Response']=='SomHit')][0] | [(S1master_log.loc[index, 'Stim/Block/Response']=='VisHit')][0]:
             S1master_log.loc[index, 'Reward']=1
        else:
            S1master_log.loc[index, 'Reward']=0
            
            
            
###############################################################################
# 'rewarded_licks','non_rewarded_licks'
###############################################################################           
        
        
#'rewarded_licks' and 'unrewarded_licks'
Name_df=['rewarded_licks','non_rewarded_licks']
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
#end of preprocessing
                

for mouse in mice:
    for i in S1master_log[S1master_log['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        if S1master_log.loc[i,'Reward']==1:
            if S1master_log.loc[i,'block_type']=='Whisker':
                a=S1master_log.loc[i,'licks_left']>0.01
                b=S1master_log.loc[i,'licks_left']<3
                mask = [all(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                S1master_log.at[i,'rewarded_licks']=S1master_log.loc[i,'licks_left'][mask] #only take the licks that are within 3secs of stim
                temp1=np.ravel(S1master_log.loc[i,'licks_left'][S1master_log.loc[i,'licks_left']<0.01])# licks from that trial before the stim are unrewarded
                temp2=np.ravel(S1master_log.loc[i,'licks_right']) # licks on wrong port
                non_rewarded_licks=np.hstack((temp1, temp2))
                S1master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
            elif S1master_log.loc[i,'block_type']=='Visual':
                a=S1master_log.loc[i,'licks_right']>0.01
                b=S1master_log.loc[i,'licks_right']<3
                mask = [all(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                S1master_log.at[i,'rewarded_licks']=S1master_log.loc[i,'licks_right'][mask] #only take the licks that are within 3secs of stim
                temp1=np.ravel(S1master_log.loc[i,'licks_right'][S1master_log.loc[i,'licks_right']<0.01])# licks from that trial before the stim are unrewarded
                temp2=np.ravel(S1master_log.loc[i,'licks_left']) # licks on wrong port
                non_rewarded_licks=np.hstack((temp1, temp2))
                S1master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
        else:
            S1master_log.at[i,'rewarded_licks']=[]
            temp1=np.ravel(S1master_log.loc[i,'licks_right'])
            temp2=np.ravel(S1master_log.loc[i,'licks_left'])
            non_rewarded_licks=np.hstack((temp1, temp2))
            S1master_log.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                            
##########################################################################################################
## Repeat the non_rewarded lick analysis but only count the licks that are not in the middle of bouts. So only isolated licks or first lick of a bout
##########################################################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licksNonBouts', 'non_rewarded_Left_licksNonBouts'    ]
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
    
for mouse in mice:
    for i in S1master_log[S1master_log['mouse_name']==mouse].index:
        temp_right=[]
        temp_left=[]
        trial_log=S1master_log.loc[i,:]
        for ind in range(0, len(trial_log.non_rewarded_licks)):
            if ind==0: #if first lick analyzed no need to check before
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right']:
                    temp_right.append(trial_log['non_rewarded_licks'][ind])
               
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left']:
                    temp_left.append(trial_log['non_rewarded_licks'][ind])
            elif trial_log['non_rewarded_licks'][ind] - trial_log['non_rewarded_licks'][ind-1]>0.5: # so only if the lick is not preceded by a lick within 500msec
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right']:
                    temp_right.append(trial_log['non_rewarded_licks'][ind])
               
                if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left']:
                    temp_left.append(trial_log['non_rewarded_licks'][ind])
        S1master_log.at[i,'non_rewarded_Right_licksNonBouts']=temp_right
        S1master_log.at[i,'non_rewarded_Left_licksNonBouts']=temp_left
        test=S1master_log['non_rewarded_Left_licksNonBouts']
        
# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
NR_Right_licks_aligned_spikes_NonBouts = []
NR_Left_licks_aligned_spikes_NonBouts = []

for num in S1master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    temp_right = []
    temp_left = []
    if np.sum(S1master_log.non_rewarded_licks[num])>0:
        for ind in range(0, len(S1master_log.non_rewarded_licks[num])):
            if S1master_log['non_rewarded_licks'][num][ind] in  S1master_log['non_rewarded_Right_licksNonBouts'][num]:
                temp_right.append(np.array([np.squeeze(S1master_log['spike_times'][num])-S1master_log['non_rewarded_licks'][num][ind]]))
            elif S1master_log['non_rewarded_licks'][num][ind] in  S1master_log['non_rewarded_Left_licksNonBouts'][num]:
                temp_left.append(np.array([np.squeeze(S1master_log['spike_times'][num])-S1master_log['non_rewarded_licks'][num][ind]]))

                
    NR_Right_licks_aligned_spikes_NonBouts.append(temp_right)
    NR_Left_licks_aligned_spikes_NonBouts.append(temp_left)
S1master_log['NR_Right_licks_aligned_spikes_NonBouts'] = NR_Right_licks_aligned_spikes_NonBouts
S1master_log['NR_Left_licks_aligned_spikes_NonBouts'] = NR_Left_licks_aligned_spikes_NonBouts
            
##########################################################################################################
## Repeat the ITBNs but adding the limitation of no lick Following" on top of no lick "preceding"
##########################################################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licksNonBouts2', 'non_rewarded_Left_licksNonBouts2'    ]
for i,name in enumerate(Name_df):   
    S1master_log[Name_df[i]]=np.zeros(len(S1master_log))
    S1master_log[Name_df[i]]=S1master_log[Name_df[i]].astype(object)
    
for mouse in mice:
    print(mouse)
    for i in S1master_log[S1master_log['mouse_name']==mouse].index:
        temp_right=[]
        temp_left=[]
        trial_log=S1master_log.loc[i,:]
        if len(trial_log.non_rewarded_licks)==1: #if there are only 1 lick, count it
            if trial_log['non_rewarded_licks'][0] in  trial_log['licks_right']:
                    temp_right.append(trial_log['non_rewarded_licks'][0])
            if trial_log['non_rewarded_licks'][0] in  trial_log['licks_left']:
                    temp_left.append(trial_log['non_rewarded_licks'][0])
        if len(trial_log.non_rewarded_licks)==2: #if there are only 2 lick, count them only if sep by 500msec
            if (trial_log['non_rewarded_licks'][1] -  trial_log['non_rewarded_licks'][0])>0.5:
                for ind in range(0, len(trial_log.non_rewarded_licks)):
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right']:
                            temp_right.append(trial_log['non_rewarded_licks'][ind])
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left']:
                            temp_left.append(trial_log['non_rewarded_licks'][ind])
        else:
            for ind in range(0, len(trial_log.non_rewarded_licks))[1:-1]:
                if ((trial_log['non_rewarded_licks'][ind] - trial_log['non_rewarded_licks'][ind-1]>0.5) & (trial_log['non_rewarded_licks'][ind+1] - trial_log['non_rewarded_licks'][ind]>0.5)): # so only if the lick isolated by 500msec
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right']:
                        temp_right.append(trial_log['non_rewarded_licks'][ind])
                       
                    if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left']:
                        temp_left.append(trial_log['non_rewarded_licks'][ind])
                    
        S1master_log.at[i,'non_rewarded_Right_licksNonBouts2']=temp_right
        S1master_log.at[i,'non_rewarded_Left_licksNonBouts2']=temp_left
        
# 'NR_Right_licks_aligned_spikes', 'NR_Left_licks_aligned_spikes'
NR_Right_licks_aligned_spikes_NonBouts = []
NR_Left_licks_aligned_spikes_NonBouts = []

for num in S1master_log.non_rewarded_licks.index:
    if num%1000 == 0:
            print(num)
    temp_right = []
    temp_left = []
    
    if np.sum(S1master_log.non_rewarded_licks[num])>0:
        for ind in range(0, len(S1master_log.non_rewarded_licks[num])):
            if S1master_log['non_rewarded_licks'][num][ind] in  S1master_log['non_rewarded_Right_licksNonBouts2'][num]:
                temp_right.append(np.array([np.squeeze(S1master_log['spike_times'][num])-S1master_log['non_rewarded_licks'][num][ind]]))
            elif S1master_log['non_rewarded_licks'][num][ind] in  S1master_log['non_rewarded_Left_licksNonBouts2'][num]:
                temp_left.append(np.array([np.squeeze(S1master_log['spike_times'][num])-S1master_log['non_rewarded_licks'][num][ind]]))

                
    NR_Right_licks_aligned_spikes_NonBouts.append(temp_right)
    NR_Left_licks_aligned_spikes_NonBouts.append(temp_left)
S1master_log['NR_Right_licks_aligned_spikes_NonBouts2'] = NR_Right_licks_aligned_spikes_NonBouts
S1master_log['NR_Left_licks_aligned_spikes_NonBouts2'] = NR_Left_licks_aligned_spikes_NonBouts
          

###############################################################################
#unit_name
###############################################################################
for mouse in mice:
    mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for cluster in np.unique(day_log['cluster_name']):
            cluster_log=day_log.loc[np.equal(day_log['cluster_name'], cluster)]
            S1master_log.loc[cluster_log.index,'unit_name']= cluster_log['mouse_name'].values[0]+'_'+ cluster_log['date'].values[0]+'_'+ cluster_log['cluster_name'].values[0]
            print(cluster_log['mouse_name'].values[0]+'_'+ cluster_log['date'].values[0]+'_'+ cluster_log['cluster_name'].values[0])


###############################################################################
#SAVE
###############################################################################
drive_path='C:\\Users\Brown Lab\Desktop\Analysis working folder\Claustrum explore'
S1master_log.to_pickle(   drive_path+'\S1master_logC.pkl', protocol=4)# for python 3


# %%
##############################################################################
##############################################################################
# Generate XCorr_df and S1_XCorr_df
##############################################################################
##############################################################################
from Claustrum_preprocessing_functions import Generate_XCorr_df
XCorr_df=Generate_XCorr_df(master_log)

from Claustrum_preprocessing_functions import Generate_S1XCorr_df
S1XCorr_df=Generate_S1XCorr_df(S1master_log)

from Claustrum_preprocessing_functions import  Add_POST_to_XCorr_df
XCorr_df= Add_POST_to_XCorr_df(master_log, XCorr_df)

from Claustrum_preprocessing_functions import  Add_POST_to_S1XCorr_df
S1XCorr_df= Add_POST_to_S1XCorr_df(S1master_log, S1XCorr_df)

#Add 'LickPref' to XCorr_df
for each in XCorr_df.index:
    Neuron1=XCorr_df.loc[each, 'Neuron1_name']
    Neuron2=XCorr_df.loc[each, 'Neuron2_name']
    Lick1=master_log[master_log['unit_name']==Neuron1]['LickPref'].values[0]
    Lick2=master_log[master_log['unit_name']==Neuron2]['LickPref'].values[0]
    XCorr_df.at[each,'Neuron1_LickPref']=Lick1
    XCorr_df.at[each,'Neuron2_LickPref']=Lick2

###############################################################################
#SAVE
###############################################################################   
XCorr_df.to_pickle('C:\\Users\Brown Lab\Desktop\master_log\FINAL\XCorr_df_20201109.pkl', protocol=4)    
S1XCorr_df.to_pickle('C:\\Users\Brown Lab\Desktop\master_log\FINAL\S1XCorr_df_20200730.pkl', protocol=4)    

# %%
##############################################################################
##############################################################################
# Generate master_DREADD
##############################################################################
##############################################################################

###############################################################################
################################## MATLAB #####################################
###############################################################################

#Run: Generate_Log_IntanData_forBehavior.m for each mouse
#Run: 'Run_Add_Trial_StartStop_to_log.m' for each mouse



# Basic initial DF construction: (Here we have only behvaior, no units)
columns=['mouse_name', 'date','block_type', 'trial_type', 'touch_stimulus',
 'vis_stimulus', 'response', 'trial_num', 'stim_onset', 'stim_offset',
 'licks_right', 'licks_left','TrialOptoStimStarts','TrialOptoStimEnds'] #these are the initial columns used to build the df, more are added as more parameters are needed
master_DREADD=pd.DataFrame(columns=columns)
#fnall=[]
#fnall += [mice for mice in os.listdir(drive_path) if mice.startswith('Claustrum')] #find all the mouse folders
#for mouse in fnall: 
#    mouse_path=os.path.join(drive_path, mouse)
drive_path = 'D:\Claustrum\DREADD'
mice=['ClaustrumO','ClaustrumP','ClaustrumQ','ClaustrumR','ClaustrumS', 'ClaustrumT', 'ClaustrumU', 'ClaustrumV','ClaustrumW','ClaustrumX', 'ClaustrumY']
for mouse in mice:
      mouse_dir=os.path.join(drive_path,mouse)
      fnall2=[]
      fnall2 += [day for day in os.listdir(mouse_dir) if day.startswith('2021')] # find all the date folders
      for date in fnall2[:8]:
          day_path=os.path.join(mouse_dir,date)
          day=[]
          day+=[matlog for matlog in os.listdir(day_path) if matlog.startswith('Log')]
          DayLog=loadmat(day_path+'/'+str(day)[2:-2])
          log=pd.DataFrame(columns=columns)
          log['mouse_name']=DayLog['log'][:,0]
          log['date']=DayLog['log'][:,1]
          log['block_type']=DayLog['log'][:,2]
          log['trial_type']=DayLog['log'][:,3]
          log['touch_stimulus']=DayLog['log'][:,4]
          log['vis_stimulus']=DayLog['log'][:,5]
          log['response']=DayLog['log'][:,6]
          log['trial_num']=DayLog['log'][:,7]
          log['stim_onset']=DayLog['log'][:,8]
          log['stim_offset']=DayLog['log'][:,9]
          log['licks_right']=DayLog['log'][:,10]
          log['licks_left']=DayLog['log'][:,11]
          log['TrialOptoStimStarts']=DayLog['log'][:,14]
          log['TrialOptoStimEnds']=DayLog['log'][:,15]
          log['Trial_Start_Stop']=DayLog['log'][:,16] #requires Add_Trial_StartStop_to_log.m to have been run
          master_DREADD=pd.concat([master_DREADD,log],ignore_index=True)
          print(mouse+date)
          
###############################################################################

#Only run the first time (initialize columns)
Name_df=['licks_right_trialonly' ,'licks_left_trialonly']
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
    
#pick which mouse you want to process
for mouse in mice:
    for i in master_DREADD[master_DREADD['mouse_name']==mouse].index:
        if i%100==0:
            print(i)
        Trial_length=master_DREADD.loc[i,'Trial_Start_Stop'][0][1] - master_DREADD.loc[i,'Trial_Start_Stop'][0][0] #get the length of trial Stop-start
        right_index=np.ravel(master_DREADD.loc[i,'licks_right'][0][0]<Trial_length) #index of right licks within trial
        left_index=np.ravel(master_DREADD.loc[i,'licks_left'][0][0]<Trial_length)#index of left licks within trial
        right_licks=master_DREADD.loc[i,'licks_right'][0][0]
        left_licks=master_DREADD.loc[i,'licks_left'][0][0]
        master_DREADD.at[i,'licks_right_trialonly']=right_licks[right_index]
        master_DREADD.at[i,'licks_left_trialonly']=left_licks[left_index]
        

        
###############################################################################
#Add 'Stim/Block/Response'
############################################################################### 
#For the Original task: Som:lick right, vis:lick left       
#For the REVERSED task: Som:lick left, vis:lick right
for mouse in mice:  
    for index in master_DREADD[(master_DREADD['mouse_name']==mouse)].index:
        if index%100==0:
            print(index) 
        if [(master_DREADD.loc[index,'block_type']=='Whisker') & (master_DREADD.loc[index, 'trial_type']=='Stim_Som_NoCue') & (master_DREADD.loc[index,'response']==1)][0]:        
            master_DREADD.loc[index, 'Stim/Block/Response']='SomHit'
        elif [(master_DREADD.loc[index,'block_type']=='Whisker') & (master_DREADD.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_DREADD.loc[index,'response']==0)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='SomMiss'
        elif [(master_DREADD.loc[index,'block_type']=='Visual') & (master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_DREADD.loc[index,'response']==2)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisHit'
        elif [(master_DREADD.loc[index,'block_type']=='Visual') & (master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_DREADD.loc[index,'response']==0)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisMiss'
        elif [(master_DREADD.loc[index,'block_type']=='Whisker') & (master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_DREADD.loc[index,'response']==0)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='SomCR'
        elif [( master_DREADD.loc[index,'block_type']=='Whisker') & ( master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_DREADD.loc[index,'response']==2)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='SomFA1'
        elif [( master_DREADD.loc[index,'block_type']=='Whisker') & ( master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & (master_DREADD.loc[index,'response']==1)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='SomFA2'
        elif [( master_DREADD.loc[index,'block_type']=='Whisker') & ( master_DREADD.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_DREADD.loc[index,'response']==2)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='SomFA3'
        elif [(master_DREADD.loc[index,'block_type']=='Visual') & (master_DREADD.loc[index,'trial_type']=='Stim_Som_NoCue') & (master_DREADD.loc[index,'response']==0)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisCR'
        elif [( master_DREADD.loc[index,'block_type']=='Visual') & ( master_DREADD.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_DREADD.loc[index,'response']==2)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisFA2'
        elif [( master_DREADD.loc[index,'block_type']=='Visual') & ( master_DREADD.loc[index,'trial_type']=='Stim_Som_NoCue') & ( master_DREADD.loc[index,'response']==1)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisFA1'
        elif [( master_DREADD.loc[index,'block_type']=='Visual') & ( master_DREADD.loc[index,'trial_type']=='Stim_Vis_NoCue') & ( master_DREADD.loc[index,'response']==1)][0]:
            master_DREADD.loc[index, 'Stim/Block/Response']='VisFA3'

################################################################################
# clean up the date brackets
###############################################################################
for each in master_DREADD.index:
    print(master_DREADD.loc[each,'date'][0])
    if len(master_DREADD.loc[each,'date'][0])==1:
        master_DREADD.at[each, 'date']=     master_DREADD.loc[each,'date'][0][0]    
        
################################################################################
#Build lick-aligned variables
################################################################################
         
#INIIALIZE COLUMN
Name_df=['LeftFirstLick', 'RightFirstLick', 'FirstLick']
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing

for mouse in mice:
    mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for trial in day_log.index:
            trial_log=day_log.loc[trial]
            if trial_log['stim_onset']: 
                if not trial_log['licks_left_trialonly'].size: #if there was no lick, make it 0 such that it'll be sorted first
                    master_DREADD.at[trial,'LeftFirstLick']=0; # I can use trial as the index directly to master_log 
                else:
                    TempLicks=trial_log['licks_left_trialonly']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                    #TempLicks=TempLicks[0][0] #get rid of extra brackets
                    if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                        master_DREADD.at[trial,'LeftFirstLick']=0;
                    else:  
                        master_DREADD.at[trial,'LeftFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count
                if not trial_log['licks_right_trialonly'].size: #if there was no lick, make it 0 such that it'll be sorted first
                    master_DREADD.at[trial,'RightFirstLick']=0; # I can use trial as the index directly to master_log 
                else:
                    TempLicks=trial_log['licks_right_trialonly']-trial_log['stim_onset']; # if there are licks, substract stim time to enable identification of first lick after stim, HOWEVER we will store then as the initial timestamp values, not as stim alinged licks
                    #TempLicks=TempLicks[0][0] #get rid of extra brackets
                    if np.max(TempLicks)<0: #if all the licks are negtive, then there was no response -> 0
                        master_DREADD.at[trial,'RightFirstLick']=0;
                    else:
                        master_DREADD.at[trial,'RightFirstLick']=np.min(TempLicks[TempLicks>0]); #Licks with negative values are before the stim, therefore don't count
#The followinf loop goes through the possibilites and populates the
for mouse in mice:
    mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for trial in day_log.index:
            trial_log=day_log.loc[trial]        
            if trial_log['stim_onset']:
                if (trial_log['LeftFirstLick']!=0) | (trial_log['RightFirstLick']!=0):
                    if (trial_log['LeftFirstLick']==0) & (trial_log['RightFirstLick']!=0):
                        master_DREADD.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                    elif (trial_log['LeftFirstLick']!=0) & (trial_log['RightFirstLick']==0):
                        master_DREADD.at[trial, 'FirstLick']=trial_log['LeftFirstLick']
                    elif trial_log['LeftFirstLick']> trial_log['RightFirstLick']:
                        master_DREADD.at[trial, 'FirstLick']=trial_log['RightFirstLick']
                    elif trial_log['LeftFirstLick']< trial_log['RightFirstLick']:
                        master_DREADD.at[trial, 'FirstLick']=trial_log['LeftFirstLick']   



###############################################################################
#'Reward'
###############################################################################       
for mouse in mice:                
    for index in master_DREADD[(master_DREADD['mouse_name']==mouse)].index:
        if index%100==0:
            print(index)
        if [(master_DREADD.loc[index,'Stim/Block/Response']=='SomHit')][0] | [(master_DREADD.loc[index, 'Stim/Block/Response']=='VisHit')][0]:
             master_DREADD.loc[index, 'Reward']=1
        else:
            master_DREADD.loc[index, 'Reward']=0
            
            


###############################################################################
# 'rewarded_licks','non_rewarded_licks'
###############################################################################           
        
#'rewarded_licks' and 'unrewarded_licks'
Name_df=['rewarded_licks','non_rewarded_licks']
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
                

for mouse in mice:
    mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]
    for day in np.unique(mouse_log['date']):
        day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
        for j,i in enumerate(day_log.index):
                if master_DREADD.loc[i,'Reward']==1:
                    if master_DREADD.loc[i,'block_type']=='Visual':
                        a=(master_DREADD.loc[i,'licks_left_trialonly'] - master_DREADD.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_DREADD.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_DREADD.loc[day_log.index[j-1],'Reward']==1:
                            b=(master_DREADD.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_DREADD.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_DREADD.at[i,'rewarded_licks']=master_DREADD.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_DREADD.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_DREADD.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_DREADD.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_DREADD.loc[i,'block_type']=='Whisker':
                        a=(master_DREADD.loc[i,'licks_right_trialonly'] - master_DREADD.loc[i,'stim_onset'])>0.01 #Here relative to stim: this should give only the First lick, since it terminated the trial in Arduino
                        if j==0:
                            b=(master_DREADD.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_DREADD.loc[day_log.index[j-1],'Reward']==1: 
                            b=(master_DREADD.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_DREADD.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_DREADD.at[i,'rewarded_licks']=master_DREADD.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_DREADD.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_DREADD.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_DREADD.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                else: #no reward
                    if master_DREADD.loc[i,'block_type']=='Visual':
                        a=(master_DREADD.loc[i,'licks_left_trialonly'])<0 # All False
                        if j==0:
                            b=(master_DREADD.loc[i,'licks_left_trialonly'])<0 # All False
                        elif master_DREADD.loc[day_log.index[j-1],'Reward']==1: 
                            b=(master_DREADD.loc[i,'licks_left_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_DREADD.loc[i,'licks_left_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_DREADD.at[i,'rewarded_licks']=master_DREADD.loc[i,'licks_left_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_DREADD.loc[i,'licks_left_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_DREADD.loc[i,'licks_right_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_DREADD.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
                    elif master_DREADD.loc[i,'block_type']=='Whisker':
                        a=(master_DREADD.loc[i,'licks_right_trialonly'])<0 # All False
                        if j==0:
                            b=(master_DREADD.loc[i,'licks_right_trialonly'])<0 # All False
                        elif master_DREADD.loc[day_log.index[j-1],'Reward']==1: 
                            b=(master_DREADD.loc[i,'licks_right_trialonly'] )<3 #here relative to trial_start: count the first 3sec, which are the 3sec following the first lick of the previous trial
                        else:
                            b=(master_DREADD.loc[i,'licks_right_trialonly'])<0 # All False
                        mask = [any(tup) for tup in zip(a, b)] # this is the right way to do boolean masks
                        master_DREADD.at[i,'rewarded_licks']=master_DREADD.loc[i,'licks_right_trialonly'][mask] #only take the licks that are within 3secs of stim
                        temp1=np.ravel(master_DREADD.loc[i,'licks_right_trialonly'][[False if x else True for x in mask]])# licks from that trial before the stim are unrewarded
                        temp2=np.ravel(master_DREADD.loc[i,'licks_left_trialonly']) # licks on wrong port
                        non_rewarded_licks=np.hstack((temp1, temp2))
                        master_DREADD.at[i,'non_rewarded_licks']=non_rewarded_licks[non_rewarded_licks!=0] #need to add this step because in the original matlab log I put a '0' when there were no licks...
              

###############################################################################
#'non_rewarded_Right_licks', 'non_rewarded_Left_licks'
############################################################################### 

#'non_rewarded_Right_licks', 'non_rewarded_Left_licks'           
Name_df=['non_rewarded_Right_licks', 'non_rewarded_Left_licks'    ]
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
    
for mouse in mice:
    for i in master_DREADD[master_DREADD['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        temp_right=[]
        temp_left=[]
        trial_log=master_DREADD.loc[i,:]
        for ind in range(0, len(trial_log.non_rewarded_licks)):
            
            if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                temp_right.append(trial_log['non_rewarded_licks'][ind])
       
            if trial_log['non_rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                temp_left.append(trial_log['non_rewarded_licks'][ind])
        master_DREADD.at[i,'non_rewarded_Right_licks']=temp_right
        master_DREADD.at[i,'non_rewarded_Left_licks']=temp_left

###############################################################################
#'rewarded_Right_licks', 'rewarded_Left_licks'  
############################################################################### 
            
#'rewarded_Right_licks', 'rewarded_Left_licks'           
Name_df=['rewarded_Right_licks', 'rewarded_Left_licks'    ]
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
for mouse in mice:
    for i in master_DREADD[master_DREADD['mouse_name']==mouse].index:
        if i%1000 == 0:
            print(i)
        temp_right=[]
        temp_left=[]
        trial_log=master_DREADD.loc[i,:]
        for ind in range(0, len(trial_log.rewarded_licks)):
            
            if trial_log['rewarded_licks'][ind] in  trial_log['licks_right'][0][0]:
                temp_right.append(trial_log['rewarded_licks'][ind])
       
            if trial_log['rewarded_licks'][ind] in  trial_log['licks_left'][0][0]:
                temp_left.append(trial_log['rewarded_licks'][ind])
        master_DREADD.at[i,'rewarded_Right_licks']=temp_right
        master_DREADD.at[i,'rewarded_Left_licks']=temp_left


###############################################################################3
## ADD A 'SALINE' 'DREADD' column
###############################################################################3
            
       
#INIIALIZE COLUMN
Name_df=['Injection']
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
            
for mouse in mice:
      mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]     
      if mouse=='ClaustrumO':
            for day in ['20210501','20210503', '20210506',  '20210509']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210502','20210504', '20210505',  '20210508']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumP':
            for day in ['20210503', '20210505',  '20210509',  '20210511']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210504', '20210506',  '20210508',  '20210510']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumQ':
            for day in ['20210504', '20210506',  '20210509',  '20210511']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210505', '20210507',  '20210508',  '20210510']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumR':
            for day in [ '20210423', '20210429','20210502', '20210504']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in [ '20210424', '20210430','20210501', '20210503']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumS':
            for day in ['20210502', '20210504', '20210507',  '20210509']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210503','20210505','20210506',  '20210508']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumT':
            for day in [ '20210512',  '20210514',  '20210518',  '20210520']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210513','20210515', '20210516']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumU':
            for day in [ '20210513',  '20210515',  '20210518', '20210520']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in [ '20210514',  '20210516',  '20210517', '20210519']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumV':
            for day in ['20210510', '20210512', '20210515','20210518']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in [ '20210511', '20210513', '20210514', '20210516']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumW':
            for day in ['20210424', '20210429','20210502', '20210513']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210425', '20210430','20210501', '20210512']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumX':
            for day in ['20210504','20210506', '20210509','20210511']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in ['20210505','20210507', '20210508','20210510']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
      if mouse=='ClaustrumY':
            for day in [ '20210511', '20210513', '20210516', '20210519']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='SALINE'
            for day in [ '20210512', '20210514', '20210517', '20210518']:
                   day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                   for trial in day_log.index:
                         master_DREADD.at[trial, 'Injection']='AGONIST21'
                   
                   
##############################################################################
## ADD CONTROLvs DREADD column
##############################################################################
                        
       
#INIIALIZE COLUMN
Name_df=['virus']
for i,name in enumerate(Name_df):   
    master_DREADD[Name_df[i]]=np.zeros(len(master_DREADD))
    master_DREADD[Name_df[i]]=master_DREADD[Name_df[i]].astype(object)
#end of preprocessing
            
for mouse in mice:
      mouse_log=master_DREADD.loc[np.equal(master_DREADD['mouse_name'], mouse)]     
      if mouse=='ClaustrumO':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'      
      if mouse=='ClaustrumP':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'  
      if mouse=='ClaustrumQ':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='CONTROL'   
      if mouse=='ClaustrumR':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='CONTROL'   
      if mouse=='ClaustrumS':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='CONTROL'   
      if mouse=='ClaustrumT':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='CONTROL'   
      if mouse=='ClaustrumU':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'   
      if mouse=='ClaustrumV':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'   
      if mouse=='ClaustrumW':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'  
      if mouse=='ClaustrumX':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='CONTROL' 
      if mouse=='ClaustrumY':
            for trial in mouse_log.index:
                  master_DREADD.at[trial, 'virus']='DREADD'   
            
            
###############################################################################
############################## SAVE
###############################################################################
            
master_DREADD.to_pickle('E:\\Claustrum\DREADD\master_DREADD.pkl')