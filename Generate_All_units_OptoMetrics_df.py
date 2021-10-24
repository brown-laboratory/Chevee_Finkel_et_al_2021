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
import pickle
sys.path.append('C:\\Users\Brown Lab\Documents\Maxime_tools')
from DEF_Get_formatted_OptoLog import Get_formatted_OptoLog

All_units_OptoMetrics_df=pd.DataFrame(columns=['frequency',
 'unit_number',
 'mouse_name',
 'date',
 'cluster_name',
 'pulses',
 'all_spikes',
 'control_spikes',
 'evoked_spikes',
 'reliability',
 'latencies',
 'mean_latencies',
 'latency_bin_starts',
 'post_pulse_zScores',
 'waveform_corr',
 'mean_cont_waveform1',
 'mean_cont_waveform2',
 'mean_cont_waveform3',
 'mean_cont_waveform4',
 'mean_evoked_waveform1',
 'mean_evoked_waveform2',
 'mean_evoked_waveform3',
 'mean_evoked_waveform4'])

mice_path=[
        'I:\Claustrum4',
        'I:\Claustrum5',
        'I:\Claustrum6']
for mouse_path in mice_path:
    os.chdir( mouse_path) # go inside the mouse folder
    fnall=[]
    fnall += [date for date in os.listdir(mouse_path) if date.endswith('-17')] #find all the date folders
    for date in fnall:
        print(date)
        date_path=os.path.join(mouse_path,date)
        os.chdir(date_path)
        Date_opto_df=Get_formatted_OptoLog()
        All_units_OptoMetrics_df=pd.concat([All_units_OptoMetrics_df,Date_opto_df])

mice_path=[
        'I:\Claustrum31',
        'I:\Claustrum32',
        'I:\Claustrum37',
        'D:\Claustrum18', 
        'D:\Claustrum23', 
        'D:\Claustrum25']
for mouse_path in mice_path:
    os.chdir( mouse_path) # go inside the mouse folder
    fnall=[]
    fnall += [date for date in os.listdir(mouse_path) if date.startswith('2019')] #find all the date folders
    for date in fnall:
        print(date)
        date_path=os.path.join(mouse_path,date)
        os.chdir(date_path)
        Date_opto_df=Get_formatted_OptoLog()
        All_units_OptoMetrics_df=pd.concat([All_units_OptoMetrics_df,Date_opto_df])
        
All_units_OptoMetrics_df=All_units_OptoMetrics_df.reset_index()
   
All_units_OptoMetrics_df['unit_name']=np.zeros(len(All_units_OptoMetrics_df))
for idx in All_units_OptoMetrics_df.index:
       All_units_OptoMetrics_df.loc[idx,'unit_name']= str(All_units_OptoMetrics_df.loc[idx,'mouse_name'])+'_'+ str(All_units_OptoMetrics_df.loc[idx,'date'])+'_'+ str(All_units_OptoMetrics_df.loc[idx,'cluster_name'])   

All_units_OptoMetrics_df.to_pickle( 'C:\\Users\Brown Lab\Desktop\master_log\All_units_OptoMetrics_df.pkl', protocol=4)
