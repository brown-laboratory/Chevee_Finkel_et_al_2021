################################################################################
# Claustrum_preprocessing_functions: contains functions called by 
# Generate_DataFrames_CheveeFinkeletal2021.py
################################################################################

import scipy as sp
import scipy.io
import scipy.stats
import os
import numpy as np
import pandas as pd
import glob
import csv
from tqdm import tnrange, tqdm_notebook
from collections import Iterable
import matplotlib.pylab as mpl
import matplotlib.patches as patches
from matplotlib import gridspec
import seaborn as sns
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}

mpl.rc('font', **font)


def opto_identify(unit_num, log_df, opto_log_df, opto_spikes_df):

    unit = opto_log_df.iloc[unit_num]
    unit_indices = (unit['cluster_name'] ==
                    opto_spikes_df['cluster_name'])
    unit_spikes = opto_spikes_df.loc[unit_indices, 'spikes']

    mpl.close('all')
    frequencies = np.array(list(unit['first_last_opto_pulses'].keys()))
    #frequencies = frequencies[frequencies != 2]

    fig = mpl.figure(figsize=(10, 20))
    fig.suptitle(unit['mouse_name'] + ', ' +  unit['date'] + ', ' + unit['cluster_name'])
    gs = gridspec.GridSpec(11, 1, height_ratios=[1, 15, 5, 1, 15, 5, 1, 15, 5, 1, 15])

    raster_positions = [1,4,7,10]
    stim_fig_positions = [0,3,6,9]

    for i, frequency in enumerate(frequencies):
        
        raster = fig.add_subplot(gs[raster_positions[i]])
        stim_fig = fig.add_subplot(gs[stim_fig_positions[i]], sharex = raster)
        
        flop = unit['first_last_opto_pulses'][frequency]

        trial_dur = flop[1][0] - flop[0][0]
        trial_total = 0
        for trial in range(len(flop[0])):
            trial_spike_inds = (flop[0][trial]-0.5 < unit_spikes) & (unit_spikes < flop[0][trial] + trial_dur + .5)
            trial_spikes = unit_spikes[trial_spike_inds] - flop[0][trial]
            raster.vlines(trial_spikes, trial + .5, trial + 1.3, linewidth = 0.5)
            trial_total += 1

        figure_pulse_inds = [(unit['grouped_opto_pulses'][frequency] >= flop[0][-1])
                  & (unit['grouped_opto_pulses'][frequency] <= flop[1][-1])]
        figure_pulses = unit['grouped_opto_pulses'][frequency][figure_pulse_inds] - flop[0][-1]
        example_pulse_ind = np.where(unit['opto_stim_onsets'] == flop[0][-1])
        stim_duration = unit['opto_stim_offsets'][example_pulse_ind] - unit['opto_stim_onsets'][example_pulse_ind]

        for p in figure_pulses:
            stim_fig.add_patch(patches.Rectangle((p,0), stim_duration, 4, color = '#3182bd')) 

        raster.autoscale(enable=True, tight=True)
        raster.spines['right'].set_visible(False)
        raster.spines['top'].set_visible(False)
        raster.xaxis.set_ticks_position('bottom')
        raster.yaxis.set_ticks_position('left')
        raster.set_xlabel('Time(s)')
        raster.set_ylabel('Trials')

        stim_fig.set_ylim(0, 4,)
        stim_fig.set_xlim(-0.5, trial_dur + 0.5)
        stim_fig.axis('off')
        stim_fig.set_title(str(frequency) + ' Hz')

    mpl.subplots_adjust(left=0.1, right=.9, top=0.9, bottom=0.1)
    fig
    return fig

def find_latencies(frequency_row, unit_spikes):
    ## defines the latencies of spikes by first defining a window around the peak of
    ## first spike histogram. Spike latencies are only considered if the spike lands
    ## within that window. Abscence of spikes within the window will degrade reliability.
    spikes = unit_spikes.as_matrix()
    pulse_bins = [np.arange(pulse, pulse+0.020, 0.001)[:20] for pulse in frequency_row['pulses']] # 1 ms bins
    hist = np.mean([np.histogram(spikes, pulse)[0] for i, pulse in enumerate(pulse_bins)], axis = 0)

    # look for spikes within +/- 3 bins (3 ms) of mean latency bin 
    mean_biggest_bins = np.convolve(hist, np.ones((3,))/3, mode='same').argsort()[-1:][0] #finds middle of largest 3 consecutive bins
    window_template = [[mean_biggest_bins - 3, mean_biggest_bins + 3] if 19-3 >= mean_biggest_bins >= 3 
                       else [0, mean_biggest_bins + 3] if mean_biggest_bins <= 19-3 
                       else [mean_biggest_bins - 3, 19]][0]
    window_starts = [pulse[window_template[0]] for pulse in pulse_bins]
    window_ends = [pulse[window_template[1]] for pulse in pulse_bins]
    latency_bin_start = window_template[0] #number of ms (1ms = 1bin) after pulse where spikes are considered 'laser evoked'

    idx = spikes.searchsorted(window_starts, side = 'right')
    idx2 = idx[idx < len(unit_spikes)] #otherwise will return an error if no spikes follow the last laser pulse
    evoked_spikes = np.array([spike if spike <= window_ends[i] else np.nan for i,spike in enumerate(spikes[idx2])])
    latencies = np.array([spike - frequency_row['pulses'][i] if spike <= window_ends[i] else np.nan
                 for i,spike in enumerate(spikes[idx2])])
    reliability = [np.count_nonzero(~np.isnan(latencies))/len(latencies) if len(latencies) > 0 else 0]
    latencies = latencies[~np.isnan(latencies)]
    evoked_spikes = evoked_spikes[~np.isnan(evoked_spikes)]
    return [[evoked_spikes], [latencies], latency_bin_start,  reliability]


def opto_metrics(uni_id, log, opto_metrics_log):

    unit = opto_metrics_log[opto_metrics_log['unit_number'] == uni_id]

    frequencies = unit['frequency'].drop_duplicates()
    stim_duration = 0.002

    mpl.close('all')
    fig2 = mpl.figure(figsize=(18, 10))
    fig2.suptitle(unit['mouse_name'].iloc[0] + ', ' +  unit['date'].iloc[0] + ', ' + unit['cluster_name'].iloc[0])

    ax1 = mpl.subplot2grid((2,2), (0,0), rowspan=2, colspan=1)
    ax2 = mpl.subplot2grid((2,2), (0,1), rowspan=1, colspan=1)
    ax3 = mpl.subplot2grid((2,2), (1,1), rowspan=1, colspan=1)

    ax4 = fig2.add_axes([0.34, 0.425, 0.1, 0.4])
    
    fig2.subplots_adjust(wspace=.4, hspace = 0.4)

    sns.despine()

    data={}
    total_stims = 0

    laser_evoked_spike_inds = []
    values = np.array([])
    labels = np.array([])
    for ind in list(unit.index):
        opto_stims = unit.loc[ind, 'pulses']
        for stim_num in range(148): #This was set to 90 which means a lot of plots I generated only show about half the actual trials... it does not affect the calling of optotag, just the plot (20190809)
            spike_inds = ((unit.loc[ind, 'all_spikes'] >= opto_stims[stim_num]-0.02) &
                         (unit.loc[ind, 'all_spikes'] < opto_stims[stim_num] + .03))
            spikes = (unit.loc[ind, 'all_spikes'][spike_inds] - opto_stims[stim_num])*1000
            ax1.vlines(spikes, total_stims + stim_num + .5, total_stims + stim_num + 1.3, linewidth = 2)
            x = total_stims+stim_num

        ax1.plot([0,stim_duration*1000], [total_stims+stim_num, total_stims+stim_num], color = 'k', linewidth = 3)
        ax1.text(-7.5, total_stims + stim_num/2, str(unit.loc[ind, 'frequency']) + ' Hz'),
        total_stims = total_stims + stim_num
        
            
        latency_ind = unit.loc[ind,'latencies']<=0.025
        values = np.append(values,unit.loc[ind,'latencies'][latency_ind])
        labels = np.append(labels,[unit.loc[ind,'frequency']] * len(unit.loc[ind,'latencies'][latency_ind]))
                            
    data_df = pd.DataFrame([values*1000, labels], index = ['Latency (ms)', 'Frequency (Hz)']).T
    
    if not len(data_df)==0: #add this if statement to prevent error when there are 0 evoked spikes, which can happen for units with very few spikes
        #added by MC 20190617
        sns.violinplot(x = 'Frequency (Hz)', y = 'Latency (ms)', data = data_df, ax=ax3, color = 'xkcd:sky blue')
    
        sns.pointplot(x = 'frequency', y = 'reliability', data = unit, ax=ax2, join = False, color = 'xkcd:sky blue')
    
        offset = 0
        waveform_corr = []
        for electrode_num in range(1,5):
            laser_evoked_waveforms_df = pd.DataFrame(np.mean(unit['mean_evoked_waveform'+str(electrode_num)]))+offset
            cont_waveforms_df = pd.DataFrame(np.mean(unit['mean_cont_waveform'+str(electrode_num)]))+offset
            ax4.plot(cont_waveforms_df, color = 'k', alpha = 0.8, linewidth = 2)
            ax4.plot(laser_evoked_waveforms_df, color = 'xkcd:sky blue',alpha = 0.8, linewidth = 3)
            offset += 0.2
        ax4.set_title('r = ' + str(unit['waveform_corr'].iloc[0]))
    
        ax1.set_xlabel('Time(ms)')
        ax1.set_ylabel('Stim trial')
        ax1.autoscale(enable=True, tight=True)
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.add_patch(patches.Rectangle((0,0), stim_duration*1000, total_stims, color = 'xkcd:sky blue'))
        ax1.set_xlim(-20,30)
    
        ax3.set_xlabel('Frequency (Hz)')
        ax3.set_ylabel('Mean spike latency (ms)')
        ax3.set_ylim(0,20)
        ax3.yaxis.set_ticks(np.arange(0, 20, 2.5))
    
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('P(spike)')
        ax2.set_ylim(0,1)
    
        ax4.patch.set_alpha(0)
        ax4.axis('off')
    
        ax5 = fig2.add_axes([0.40, 0.40, 0.05, 0.1])
        ax5.spines['left'].set_visible(False)
        ax5.spines['top'].set_visible(False)
        ax5.yaxis.set_ticks_position('right')
        ax5.spines['right'].set_smart_bounds(True)
        ax5.spines['bottom'].set_smart_bounds(True)
        mpl.yticks((0,0.2,0.2))
        mpl.xticks((0,0.5,0.5))
        ax5.set_ylim(0,0.2)
        ax5.set_xlim(0,0.5)
        ax5.patch.set_alpha(0)
        ax5.set_xlabel('ms')
        ax5.set_ylabel('mV')
        ax5.yaxis.set_label_position('right')

    return fig2


def Wrap_function_optoID():
      from Claustrum_preprocessing_functions import opto_identify
      from Claustrum_preprocessing_functions import find_latencies
      from Claustrum_preprocessing_functions import opto_metrics
      task_data_files = glob.glob("Log_*")
      opto_data_files = glob.glob("Opto_log_*")
      opto_spike_files = glob.glob("optoSpikes_log_*")
      opto_wave_files = glob.glob("waveform_log_*")
      
      column_names1 =['mouse_name', 'date','block_type', 'trial_type', 'touch_stimulus',
                       'vis_stimulus', 'response', 'trial_num', 'stim_onset', 'stim_offset', 
                       'licks_right', 'licks_left', 'spike_times', 'cluster_name' ]
      column_names2 = ['mouse_name', 'date', 'cluster_name', 'opto_stim_onsets','opto_stim_offsets']
      column_names3 = np.concatenate((['mouse_name', 'date', 'cluster_name', 'spikes'],
                                       ['waveform_'+str(i) for i in range(128)]))
      
      log_df = pd.DataFrame([], columns = column_names1)
      opto_log_df = pd.DataFrame([], columns = column_names2)
      opto_spikes_df = pd.DataFrame()
      opto_waves_df = pd.DataFrame()
      
      for file_num in tnrange(len(glob.glob("Log_*"))):
           mat = sp.io.loadmat(task_data_files[file_num])
           mat2 = sp.io.loadmat(opto_data_files[file_num])
          
           log = mat['log']
           log2 = mat2['optoTable']
      
           indv_log_df = pd.DataFrame(log, columns = column_names1)
      
           log_df = pd.concat([log_df,indv_log_df])
           opto_log_df = pd.concat([opto_log_df,pd.DataFrame(log2, columns = column_names2)])
      
      #MAY NEED TO CHANGE BETWEEN 0 AND 1, CHECK WHICH IS THE .csv FILE
      opto_waves_df = pd.concat([opto_waves_df, pd.read_csv(opto_wave_files[1])])
      opto_spikes_df = pd.concat([opto_spikes_df, pd.read_csv(opto_spike_files[0])])
      
      opto_log_df = opto_log_df.reset_index(drop = True) 
      opto_spikes_df = opto_spikes_df.reset_index(drop = True)
      
       ##spikes/waves were listed in opto_waves_df by finding all spikes within 0.025s of laser pulse
       ##this can cause spikes to be listed more than once when frequency of laser stim was >50hz
       ##This will remove double listed spikes:
      opto_waves_df = opto_waves_df.drop_duplicates().reset_index(drop = True)
                           
       ##......................................
       #unique_sessions = opto_log_df[['mouse_name', 'date']].drop_duplicates().reset_index(drop = True)
      
      opto_log_df['first_last_opto_pulses'] = np.nan
      opto_log_df['grouped_opto_pulses'] = np.nan
      
      #for session in tnrange(unique_sessions.shape[0]):
      session=0 #because I will run it independently for each session
      session_row_ind = np.arange(len(opto_log_df))
      rows = opto_log_df
          
      ISIs = np.around(rows.iloc[0,3][0][1:]-rows.iloc[0,3][0][0:-1], 4)
      
      ISIs = np.concatenate(([ISIs[0]], ISIs, [ISIs[-1]]))
      unique_ISIs = np.unique(ISIs, return_counts = True)
      ind = unique_ISIs[1]>50
      unique_ISIs = unique_ISIs[0][ind]
      
      first_last_opto_pulses = {}
      grouped_opto_pulses = {}
      for isi in unique_ISIs:
        ISI_category = np.where(np.absolute(ISIs-isi) <= 0.001)[0]
        last_ind = np.where(np.diff(ISI_category) > 30)
        last_ind = last_ind[0]
        if last_ind.size ==0:
            last_ind = ISI_category[-1]
        opto_pulse_inds= np.array(range(np.min(ISI_category),last_ind))
        grouped_opto_pulses[1/isi] = rows.iloc[0,3][0][opto_pulse_inds]
      
        IBI_ind = ISIs[opto_pulse_inds] != isi
        first_pulse_ind = np.concatenate(([opto_pulse_inds[0]],opto_pulse_inds[IBI_ind]))
        last_pulse_ind = np.concatenate((first_pulse_ind[1:]-1, [opto_pulse_inds[-1]]))
        first_last_opto_pulses[1/isi] = [rows.iloc[0,3][0][first_pulse_ind],rows.iloc[0,3][0][last_pulse_ind]]
      opto_log_df.loc[session_row_ind,'first_last_opto_pulses'] = [first_last_opto_pulses]
      opto_log_df.loc[session_row_ind,'grouped_opto_pulses'] = [grouped_opto_pulses]
       ##..................................
      #  for col in [0,1,2,3,7,7,7,7,8,8,8,8,10,10,11,11,13]:
      #         log_df.iloc[:,col] = log_df.iloc[:,col].str[0]
      for col in [0,1,2,3,4]:
               opto_log_df.iloc[:,col] = opto_log_df.iloc[:,col].str[0]
      opto_spikes_df['cluster_name'] = opto_spikes_df['cluster_name'].apply(lambda y: 'T'+ y)
      opto_waves_df['cluster_name'] = opto_waves_df['cluster_name'].apply(lambda y: 'T'+ y)
       
      
      opto_log_df.head()
      
      opto_waves_df.head()


      try:
          all_opto_metrics_df = pd.read_hdf('all_opto_metrics_df.h5', 'table')
      except:
          all_opto_metrics_df = pd.DataFrame()
          
      ##.............................................................................
      ##........................................................................
          #missing_sessions = []
      for row in tnrange(opto_log_df.shape[0]):
          unit = opto_log_df.iloc[row,:]
          
      #    if all_opto_metrics_df.shape[0] != 0:
      #        if (unit[['mouse_name','date','cluster_name']] == 
      #            all_opto_metrics_df[['mouse_name','date','cluster_name']]).all(axis = 1).any():
      #            continue
                  
          df1 = unit['cluster_name']
          df2 = opto_spikes_df['cluster_name']
          df3 =  opto_waves_df['cluster_name']
          
          unit_spike_indices = pd.eval('df1 == df2')
          unit_wave_indices = pd.eval('df1 == df3')
          
      #    if not any(unit_wave_indices):
      #        missing_sessions.append([unit[['mouse_name', 'date']].as_matrix()])
      #        continue
      
          unit_spikes = opto_spikes_df.loc[unit_spike_indices, 'spikes'].reset_index(drop = True)
          unit_waves = opto_waves_df.loc[unit_wave_indices, :].reset_index(drop = True)
      
          ## create a small df to encapsulate all the data for this unit, fill in all the
          ## identifying information for the unit. Create rows for data pertaining to each
          ## optogenetic stimulus frequency given
          
          if len(unit['opto_stim_onsets']) == 0:
              
              frequencies = np.array([1,5,10,40])   
              opto_metrics_df= pd.DataFrame(frequencies, columns = ['frequency'])
              opto_metrics_df['unit_number'] = row
              opto_metrics_df['mouse_name'] = unit['mouse_name']
              opto_metrics_df['date'] = unit['date']
              opto_metrics_df['cluster_name'] = unit['cluster_name']
              opto_metrics_df['evoked_spikes'] = np.nan
              opto_metrics_df['reliability'] = np.nan
              opto_metrics_df['latencies'] = np.nan
              opto_metrics_df['mean_latencies'] = np.nan
              opto_metrics_df['latency_bin_starts']= np.nan
              opto_metrics_df['all_spikes'] = np.nan
              opto_metrics_df['control_spikes'] = np.nan
              opto_metrics_df['waveform_corr'] = np.nan
              opto_metrics_df['mean_cont_waveform1'] = np.nan
              opto_metrics_df['mean_cont_waveform2'] = np.nan
              opto_metrics_df['mean_cont_waveform3'] = np.nan
              opto_metrics_df['mean_cont_waveform4'] = np.nan
              opto_metrics_df['mean_evoked_waveform1'] = np.nan
              opto_metrics_df['mean_evoked_waveform2'] = np.nan
              opto_metrics_df['mean_evoked_waveform3'] = np.nan
              opto_metrics_df['mean_evoked_waveform4'] = np.nan
          
          frequencies = np.array(list(unit['grouped_opto_pulses'].keys()))
          #frequencies = frequencies[frequencies != 2]
      
          opto_metrics_df= pd.DataFrame(frequencies, columns = ['frequency'])
          opto_metrics_df['unit_number'] = row
          opto_metrics_df['mouse_name'] = unit['mouse_name']
          opto_metrics_df['date'] = unit['date']
          opto_metrics_df['cluster_name'] = unit['cluster_name']
          #opto_metrics_df['uni_id'] = unit['uni_id']
      
          ## select the times of all the optogenetic stimuli
          first_pulse = np.min(np.concatenate(list(unit['grouped_opto_pulses'].values())))
          last_pulse = np.max(np.concatenate(list(unit['grouped_opto_pulses'].values())))
          pulses = pd.DataFrame([[unit['grouped_opto_pulses'][freq]] for freq in frequencies],
                        columns = ['pulses'])
          opto_metrics_df['pulses'] = pulses
      
          ## create columns listing all spikes that occured around the time of the optogenetic
          ## stimuli presentations (for plotting rasters), 
          # %timeit evoked_spike_df = opto_windows_df.applymap(lambda y: unit_spikes[(y <= unit_spikes) & (y+0.01>= unit_spikes)].as_matrix())
      
          all_spikes = {}
          control_spikes = {}
          for i, freq in enumerate(frequencies):
              control_spikes[i] = unit_waves['spikes'][unit_waves['spikes'] <= first_pulse]
              all_spikes[i] = np.array(unit_spikes[(opto_metrics_df.loc[i, 'pulses'][0]-1 < unit_spikes) &
                                                   (opto_metrics_df.loc[i, 'pulses'][-1]+1 > unit_spikes)])
              opto_metrics_df['all_spikes'] = pd.Series(all_spikes)
              opto_metrics_df['control_spikes'] = pd.Series(all_spikes)
      
      
          ## list the latencies of the first spike that occured following each laser stimulus 
          latency_metrics = opto_metrics_df.apply(lambda y: find_latencies(y,unit_spikes), axis = 1) #find_latencies function defined above
          opto_metrics_df['evoked_spikes'] = [latency_metrics[i][0][0] for i in range(len(latency_metrics))]
          opto_metrics_df['reliability'] = [latency_metrics[i][3][0] for i in range(len(latency_metrics))]
          opto_metrics_df['latencies'] = [latency_metrics[i][1][0] for i in range(len(latency_metrics))]
          opto_metrics_df['mean_latencies'] = opto_metrics_df['latencies'].apply(lambda y: np.mean(y) if len(y)>0 else None)
          opto_metrics_df['latency_bin_starts']= [latency_metrics[i][2] for i in range(len(latency_metrics))]
      
          #calculate FR mean and std for Z-scores
          #all sessions have 10Hz stim- will use 1 sec before first pulse of each 10hz laser train as baseline period
          post_pulse_z_scores = {}
          for i, freq in enumerate(frequencies):
              bin_size = 0.010
              baseline_window_size = [1.5 if freq != 5 else 0.5]
              flop = unit['first_last_opto_pulses'][freq]
              baseline_bin_starts = flop[0]-baseline_window_size
              baseline_psth = np.mean([np.histogram(unit_spikes, np.arange(start,start+1.5, bin_size))[0]
                                       for start in baseline_bin_starts], axis = 0)/bin_size
              baseline_mean_std = [np.mean(baseline_psth), np.std(baseline_psth)]
              all_pulses_in_freq = unit['grouped_opto_pulses'][freq]
              latency_bin_starts = opto_metrics_df['latency_bin_starts'][i]
              test_bin_starts = [all_pulses_in_freq[(all_pulses_in_freq  >= flop[0][trial]) 
                                                   & (all_pulses_in_freq  <= flop[1][trial])] + latency_bin_starts
                                 for trial in range(len(flop[0]))]
              test_bin_starts = pd.DataFrame([window[0:10] for window in test_bin_starts]) #just looking at first 10 pulses of each trial
              ##6ms window defined in find_latency():
              test_bin_spike_counts = test_bin_starts.applymap(lambda y:
                                                               np.count_nonzero((unit_spikes >= y) 
                                                                                & (unit_spikes < y+0.006))) 
              post_pulse_z_scores[i] = (((test_bin_spike_counts.mean(axis = 0)/bin_size) 
                                         - baseline_mean_std[0])/baseline_mean_std[1]).as_matrix()
          opto_metrics_df['post_pulse_zScores'] = pd.Series(post_pulse_z_scores)
      
          all_evoked_spikes = np.concatenate(opto_metrics_df['evoked_spikes'].as_matrix())
      
          f2 = (lambda y: np.any(np.absolute(all_evoked_spikes-y) < 0.001))
          f2v = np.vectorize(f2)
          evoked_wave_ind = np.where(f2v(unit_waves['spikes']))[0]
          control_wave_ind = np.where(unit_waves['spikes'] <= first_pulse)[0]
      
          avg_control_wave = unit_waves.iloc[control_wave_ind, 4:-1].mean(axis = 0)
          avg_evoked_wave = unit_waves.iloc[evoked_wave_ind, 4:-1].mean(axis = 0)
      
          waveform_corr = [sp.stats.pearsonr(avg_control_wave, avg_evoked_wave)][0][0]
          opto_metrics_df['waveform_corr'] = np.around(pd.Series([waveform_corr]*4),2)
      
          opto_metrics_df['mean_cont_waveform1'] = pd.Series([avg_control_wave[0:32].as_matrix()]*4)
          opto_metrics_df['mean_cont_waveform2'] = pd.Series([avg_control_wave[32:64].as_matrix()]*4)
          opto_metrics_df['mean_cont_waveform3'] = pd.Series([avg_control_wave[64:96].as_matrix()]*4)
          opto_metrics_df['mean_cont_waveform4'] = pd.Series([avg_control_wave[96:].as_matrix()]*4)
          opto_metrics_df['mean_evoked_waveform1'] = pd.Series([avg_evoked_wave[0:32].as_matrix()]*4)
          opto_metrics_df['mean_evoked_waveform2'] = pd.Series([avg_evoked_wave[32:64].as_matrix()]*4)
          opto_metrics_df['mean_evoked_waveform3'] = pd.Series([avg_evoked_wave[64:96].as_matrix()]*4)
          opto_metrics_df['mean_evoked_waveform4'] = pd.Series([avg_evoked_wave[96:].as_matrix()]*4)
      
      
      
      
          # Finally concatenate to the larger dataframe
          if all_opto_metrics_df is None: 
              all_opto_metrics_df = opto_metrics_df
          else:
              all_opto_metrics_df = all_opto_metrics_df.append(opto_metrics_df)
          all_opto_metrics_df = all_opto_metrics_df.reset_index(drop = True)
      ##.............................................................................
      ##.............................................................................
          all_opto_metrics_df.head()
      ## RUN EVERYTHING BEFORE THIS IN ONE BLOCK
      frequencies = [10]
      group_names = ['identified', 'unidentified']
      
      latency_limit = 0.007
      waveform_corr_limit = 0.9
      reliability_lim = 0.3
      
      for i, freq in enumerate(frequencies):
          rows = all_opto_metrics_df['frequency'] == freq
          identified_group = ((all_opto_metrics_df['mean_latencies'][rows] <= latency_limit) &
                              (all_opto_metrics_df['waveform_corr'][rows] >= waveform_corr_limit) &
                              (all_opto_metrics_df['reliability'][rows] >= reliability_lim))
          groups = [identified_group, ~identified_group]
      
      
      
      all_units_df = all_opto_metrics_df[rows]
      identified_df = all_units_df[identified_group.as_matrix()]
      non_identified_df = all_units_df[~identified_group.as_matrix()]
      # non_identified_df = 
      non_identified_df
      
      
      from matplotlib.backends.backend_pdf import PdfPages
      mpl.close('all')
      font = {'family' : 'Arial',
      'weight' : 'normal',
      'size'   : 18}
      
      mpl.rc('font', **font)
      with PdfPages('optotest.pdf') as pdf:
          for uni_id in all_units_df.loc[:,'unit_number']:
              fig = opto_metrics(uni_id, opto_log_df, all_opto_metrics_df)
              pdf.savefig(fig)
    
      #Create proper names to be saved as Optotagged, SameTT. etc.
      FullNamesOpto = open("OptoTagged_list.txt","a") #create the text file where we will write the list of optotagged neurons, append mode
      optoTT_list=[] #create a variable where we will store the name/number of each tetrode that has a tagged neuron, we will use it to get the SameTT
      for row in identified_df.iterrows():
          Mouse=identified_df.loc[row[0],'mouse_name']
          Date=identified_df.loc[row[0],'date']
          TT=identified_df.loc[row[0],'cluster_name']
          FullNamesOpto.write(Mouse+'_'+Date+'_'+TT+'\n')
          optoTT_list.append(TT[2])
      FullNamesOpto.close()
    
      FullNamesSameTT = open("SameTT_list.txt","a") #create the text file where we will write the list of SameTT neurons, append mode
      FullNamesRest = open("Rest_list.txt","a") #create the text file where we will write the list of Rest neurons (not optotagged, not SameTT), append mode
      for row in non_identified_df.iterrows():
          Mouse=non_identified_df.loc[row[0],'mouse_name']
          Date=non_identified_df.loc[row[0],'date']
          TT=non_identified_df.loc[row[0],'cluster_name']
          if TT[2] in optoTT_list:
              FullNamesSameTT.write(Mouse+'_'+Date+'_'+TT+'\n')
          else:
              FullNamesRest.write(Mouse+'_'+Date+'_'+TT+'\n')
      FullNamesSameTT.close()
      FullNamesRest.close()

      return
  
###############################################################################
###############################################################################
######################### AUC #################################################
###############################################################################
###############################################################################
      
def Add_AUC_to_master_log(master_log):
    import sys
    mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37', 'Claustrum18', 'Claustrum23', 'Claustrum25'] 
    sys.path.insert(0,'C:\\Users\Brown Lab\Documents\Maxime_tools')
    from Helper_functions import Resample_list
    from sklearn import metrics

    
    Trial_types=[np.array(['SomHit','VisCR']),
                 np.array(['VisHit','SomCR'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    unit_SomHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset
                            
    ###############################################################################
    ###############################################################################
    # Same for block type
    ###############################################################################
    
    Trial_types=[np.array(['Whisker','Visual'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock']
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock']
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    unit_SomHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['block_type']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['block_type']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock']
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock']
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset
                            
    ###############################################################################
    ###############################################################################
    
    ###############################################################################
    # Same for Lick direction (following DHO comments 2020-10-07)
    ###############################################################################
    
    Trial_types=[np.array(['SomHit','VisHit'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit']
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit']
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    unit_SomHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit']
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit']
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset
                            
    ###############################################################################
    ###############################################################################
    ###############################################################################
    # Same for Lick direction (following DHO comments 2020-10-07)
    # LICK ALIGNED
    ###############################################################################
    
    Trial_types=[np.array(['SomHit','VisHit'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit_LickAligned']
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit_LickAligned']
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    unit_SomHit=unit_log['-2to2sec_25msecbins_LickAligned'][ (unit_log['Stim/Block/Response']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-2to2sec_25msecbins_LickAligned'][ (unit_log['Stim/Block/Response']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit_LickAligned']
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomHit_vs_VisHit_LickAligned']
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset                        
                            
    ###############################################################################
    # Same for block type, but remove trials with licks in 1sec prestim
    ###############################################################################
    
    Trial_types=[np.array(['Whisker','Visual'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock_NoLickPre']
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock_NoLickPre']
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    #REMOVE TRIALS WITH PRESTIM LICKS
                    for trial in unit_log.index:
                        trial_log=unit_log.loc[trial,:]
                        maskR=trial_log['licks_right_trialonly']-trial_log['stim_onset']
                        maskR=[1 if ((x<0) & (x>-1)) else 0 for x in maskR]
                        maskL=trial_log['licks_left_trialonly']-trial_log['stim_onset']
                        maskL=[1 if ((x<0) & (x>-1)) else 0 for x in maskL]
                        if sum(maskR+maskL):
                            unit_log=unit_log.drop(trial)
                                        
                    unit_SomHit=unit_log['-1to3sec_25msecbins_StimAligned'][(unit_log['block_type']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['block_type']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        # Because I dropped the rows of each unit_log that was to not be counte, there are values missing for these trials
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock_NoLickPre']
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_SomBlock_vs_VisBlock_NoLickPre']
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']): 
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    #need to drop the trials that we didn't ount earlier
                    index=unit_log.index
                    mask=unit_log[Column_name_95CI]!=0
                    unit_log=unit_log.loc[[x for x,y in zip(index.values,mask.values ) if y],:]
                    #
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                #print(First_Onset)
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        # Need to redifne unit_log as the full set of trials, including the ones that were dropped
                        unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                        df_index=unit_log.index
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset
                            
    ###############################################################################
    ###############################################################################     
    #Response to reviewers 2020-12-31)
    Trial_types=[np.array(['SomCR','VisCR'])]
    
    for Trial_type in Trial_types:
        #preprocessing to store the AUC values in S1master_log...
        Column_name_mean=['unit_mean_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        Column_name_95CI=['unit_95CI_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        
        for i,name in enumerate(Column_name_mean):
            
            master_log[Column_name_mean[i]]=np.zeros(len(master_log))
            master_log[Column_name_95CI[i]]=np.zeros(len(master_log))
            master_log[Column_name_mean[i]]=master_log[Column_name_mean[i]].astype(object)
            master_log[Column_name_95CI[i]]=master_log[Column_name_95CI[i]].astype(object)
        ##end of preprocessing (creating columns of the right type
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    print(unit_log['unit_name'].values[0])
        #            #now separate the trials from that unit based on what you want to decode
                    
                    unit_SomHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[0]) ]
                    X_SomHit = [i for i in unit_SomHit]
                    unit_VisHit=unit_log['-1to3sec_25msecbins_StimAligned'][ (unit_log['Stim/Block/Response']==Trial_type[1])]
                    X_VisHit = [i for i in unit_VisHit]
                    
                    
                     #create the pairs of comparisons you want to make
                    Comparisons=[(np.asarray(X_SomHit),np.asarray( X_VisHit))]
                 
                    #Need to add check that no categoris are empty beause it'll create an error (and there is no point predicting)
                    if not all([X_SomHit,  X_VisHit]):
                        print('One empty category, skip ' + np.unique(unit_log['unit_name']))
                        for i in unit_log.index:
                            for column_name_idx,Comp in enumerate(Comparisons):
                                master_log.at[i, Column_name_mean[column_name_idx]]=np.zeros((1,159))[0]
                                master_log.at[i, Column_name_95CI[column_name_idx]]=[np.zeros((1,159))[0], np.zeros((1,159))[0]]
                        continue
                        
                        
                   
                     #Loop through each category group you put in th list to be decoded
                    for column_name_idx,Comp in enumerate(Comparisons):
                        X=np.vstack(Comp)
                        y=[]
                        for i in np.arange(len(Comp[0])): #number of trials in first category
                            y.append(True)
                        for i in np.arange(len(Comp[1])): #number of trials in second category
                            y.append(False)
        
        
                        unit_array=X # extract array with spike count bins
                        unit_licks=y # extract array with the response of each trial
                        XY=np.hstack((X,np.asarray(y).reshape(len(y),1)))
                        X_reps=[]
                        for j in np.arange(1000): #Repeat the shuffle 1000 times
                            unit_array=Resample_list(XY)
                            unit_array=np.asarray(unit_array)    
                            unit_arrayX=unit_array[:,0:159]
                            unit_arrayY=unit_array[:,159].reshape(len(y),1)
                            
                            #np.random.shuffle(unit_licks);#THIS IS THE ONLY LINE THAT IS DIFFERENT THAN THE ACTUAL ROC ANALYSIS! it shuffles "in place", different random everytime
                            unit_bin_scores=np.zeros_like(unit_arrayX[0], dtype='float') #create a vector of size #bins to store the score of each bin (gathered from across trials below)
                        #Now decode each bin sequentially and store the results
                            for each in np.arange(len(unit_arrayX[0])):    #for each bin
                                scores = unit_arrayX[:,each].reshape((len(unit_arrayX),1))# get the proper bin from all trials, also some formatting...
                                auc = metrics.roc_auc_score(unit_arrayY, scores) # get AUC scores from ROC based on 
                                unit_bin_scores[each]=float(auc)
                            X_reps.append(unit_bin_scores) #at the end of each repetation, you add the result (a vector of length #bins) - ie. the predictive value of each bin with the shuffle labels.
                        X_mean=np.mean(X_reps, axis=0) # this averages the score of each bin across all 500 repetitions
                        X_95CI=[np.sort(X_reps, axis=0)[24,:], np.sort(X_reps, axis=0)[974,:]]    
                        print(Column_name_mean[column_name_idx])
                        
                        for i in unit_log.index:
                            master_log.at[i, Column_name_mean[column_name_idx]]=X_mean
                            master_log.at[i, Column_name_95CI[column_name_idx]]=X_95CI
        
        
        
        
                        
           
        Column_name_Sig=['unit_Sig_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        Column_name_Onset=['unit_Onset_25msecAUC_-1to3sec_stimAligned_'+Trial_type[0]+'_vs_'+Trial_type[1]]
        #
        for i,name in enumerate(Column_name_Sig):
            
            master_log[Column_name_Sig[i]]=np.zeros(len(master_log))
            master_log[Column_name_Onset[i]]=np.zeros(len(master_log))
            master_log[Column_name_Sig[i]]=master_log[Column_name_Sig[i]].astype(object)
            master_log[Column_name_Onset[i]]=master_log[Column_name_Onset[i]].astype(object)
        
        for mouse in mice:
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            for day in np.unique(mouse_log['date']):
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                for unit in np.unique(day_log['cluster_name']):
                    unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                    df_index=unit_log.index
                    for j,comp in enumerate(Column_name_95CI):
                        #comp is each comparison, we take the 95CI>0.5 as significant. The values in Column_name_95CI are 2 arrays, the first is top 95CI, the seconf is bottom
                        Onset_POS=unit_log.loc[df_index[0], comp][1] < 0.5 #gets a boolean vector where each bin that is significant is a one (positive scores) note: you can take the first trial, they are all the same since these values are unit specific
                        Onset_NEG=unit_log.loc[df_index[0], comp][0] > 0.5 #gets a boolean vector where each bin that is significant is a one (negative scores) note: you can take the first trial, they are all the same since these values are unit specific
                        All_Onset=Onset_POS+Onset_NEG
                        for idx, each in enumerate(All_Onset[:-3]):
                            if each & All_Onset[idx+1] & All_Onset[idx+2]: # ie. if 3 Trues back to back
                                First_Onset=np.asarray(idx) #get the index of the first bin that was significant]
                                Significant=1
                                break
                            else: # this wat, each non 3consec makes significant=0. If TTT is found, significant =1 and break
                                First_Onset=np.empty(0)
                                Significant=0
                        
                        for idx in df_index:    #the value is for the unit, but you need to populated each trial that belongs to that unit
                            master_log.at[idx,Column_name_Sig[j]]=Significant
                            master_log.at[idx,Column_name_Onset[j]]= First_Onset

    return master_log

###############################################################################
###############################################################################
######################### Response clustering #################################
###############################################################################
###############################################################################
    
def clustering_add_test(master_log):
    mice=['Cl4','Cl5','Cl6','Claustrum31', 'Claustrum32', 'Claustrum37', 'Claustrum18', 'Claustrum23', 'Claustrum25'] 
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
    Cluster_assignment=hierarchy.fcluster(Z, 5, criterion='distance')
    
    Name_df=['test']

    for i,name in enumerate(Name_df):   
        master_log[Name_df[i]]=np.zeros(len(master_log))
        master_log[Name_df[i]]=master_log[Name_df[i]].astype(object)
    #end of preprocessing
    
    for mouse in mice:
        mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
        for day in np.unique(mouse_log['date']):
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            for unit in np.unique(day_log['cluster_name']):
                print(unit)
                unit_log=day_log.loc[np.equal(day_log['cluster_name'], unit)]
                
                for j,each in enumerate(Correlations['unit_name'].index): #have to be careful with indexing
                    if Correlations.loc[each, 'unit_name'] == np.unique(unit_log['unit_name']):# find which row of Correlations corresponds to the unit
                        Cluster=Cluster_assignment[j]# get the cluster assignment
                        break #you can break out of this loop once you find the match
                for i in unit_log.index:
                    master_log.at[i, 'test']=Cluster #insert the cluster assignment for each row of master_log that corresponds to that unit
                    

    return master_log

###############################################################################
###############################################################################
######################### XCorr #################################
###############################################################################
###############################################################################

def Generate_XCorr_df(master_log)
    ORIGINAL= ['Cl4','Cl5','Cl6','Claustrum31','Claustrum32','Claustrum37']
    REVERSED= ['Claustrum18','Claustrum23','Claustrum25']
    
    mice=['Cl4','Cl5','Cl6','Claustrum18','Claustrum23','Claustrum25','Claustrum31','Claustrum32','Claustrum37']
    
    XCorr_df=pd.DataFrame(columns=['Neuron1_name', 'Neuron2_name','Neuron1_LickPref','Neuron2_LickPref', 'Neuron1_Category', 'Neuron2_Category', 'Neuron1_group', 'Neuron2_group', 'XCorr-50to50_ALL1sec', 'XCorr-50to50_shuffle_ALL1sec','XCorr-50to50_RIGHT1sec', 'XCorr-50to50_shuffle_RIGHT1sec','XCorr-50to50_LEFT1sec', 'XCorr-50to50_shuffle_LEFT1sec', 'XCorr-50to50_CORRECT1sec', 'XCorr-50to50_shuffle_CORRECT1sec', 'XCorr-50to50_INCORRECT1sec', 'XCorr-50to50_shuffle_INCORRECT1sec'])
    
    for trial_set in ['All','Right','Left','Correct', 'Incorrect']:
        Counting_pairs=0
        master_counter=0
        master_LickPref=[]
        for mouse in mice:
            print(mouse)
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            day_count=0
            for day in np.unique(mouse_log['date']):   
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                #Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num']) #these are a mess of either ints or ints trapped in arrays...
    
                if trial_set=='All':
                        Trial_numbers=np.unique(day_log['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='Correct':
                        Trial_numbers=np.unique(day_log[(day_log['Stim/Block/Response']=='SomCR') | (day_log['Stim/Block/Response']=='VisCR') | (day_log['Stim/Block/Response']=='SomHit') | (day_log['Stim/Block/Response']=='VisHit')]['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='Incorrect':
                        Trial_numbers=np.unique(day_log[(day_log['Stim/Block/Response']=='SomMiss') | (day_log['Stim/Block/Response']=='VisMiss') | (day_log['Stim/Block/Response']=='SomFA') | (day_log['Stim/Block/Response']=='VisFA')]['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='Right':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['block_type']=='Visual']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='Left':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['block_type']=='Visual']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)        
                Units_number=len(np.unique(day_log['unit_name']))
                        
                Units_number=len(np.unique(day_log['unit_name']))
                
                Lick_prefs=[]
                Names_for_plot=[]
                Category_for_plot=[]
                #First loop through and get their LickPref
                for unit in np.unique(day_log['unit_name']):
                    unit_log=day_log[day_log['unit_name']==unit]
                    Names_for_plot.append(unit_log['unit_name'].values[0])
                    Lick_prefs.append( unit_log['LickPref'].values[0])
    
                    if unit_log['Category'].values[0]=='OptoTag':
                        Category_for_plot.append('OptoTag')
                    elif unit_log['Category2'].values[0]=='OptoNetwork':
                        Category_for_plot.append('SALT')
                    elif unit_log['Category'].values[0]=='SameTT':
                        Category_for_plot.append('SameTT')
                    elif unit_log['Category3'].values[0]=='OptoNetwork_SameTT':
                        Category_for_plot.append('SALT_SameTT')
                    elif unit_log['Category4'].values[0]=='InBetween':
                        Category_for_plot.append('InBetween')
                    else:
                        Category_for_plot.append('X')
                
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                day_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in day_df.index:
                    for col in day_df.index:
                        day_df.at[row,col]=np.zeros((1,1000))
                        
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                shuffle_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in shuffle_df.index:
                    for col in shuffle_df.index:
                        shuffle_df.at[row,col]=np.zeros((1,1000))
                        
                Number_of_counted_trials=   len(Trial_numbers)-1 #we will discount each trial we remove to check that we have enough and then normalize by the right amount     
                temp=np.repeat(Number_of_counted_trials,Units_number)
                Number_of_counted_trials_PAIRMATRIX=np.asarray([ temp for _ in np.arange(Units_number)])
                for m,trial in enumerate(Trial_numbers[:-1]): #here we need to enumerate to grab the next trial when doing the control (instead of shuffle) (go until the second to last because we are taking '+1')
                    Units_log=day_log[day_log['trial_num']==trial] #represent one trial, on row per neuron
                    Units_log_NextTrial=day_log[day_log['trial_num']==Trial_numbers[m+1]] #will be used for shifted XCorr
                    
                    #Deal with error when a trial has no 'Stim_onset'               
                    a=np.unique(Units_log['stim_onset'])
                    if np.sum(a)==0:
                        a=0
                    else:
                        a=a[0][0][0][0][0]
                    b=np.unique( Units_log_NextTrial['stim_onset'])
                    if np.sum(b)==0:
                        b=0
                    else:
                        b=b[0][0][0][0][0]
    
                    if ((a<1) | (b<1)):
                        Number_of_counted_trials_PAIRMATRIX[:,:]-=1
                        print(a)
                        continue
                    
                    
                    
                    for i,neuron1 in enumerate(Units_log.index):
                        spikes1=Units_log.loc[neuron1,'spike_times']
                        spikes1=spikes1[spikes1<Units_log.loc[neuron1,'stim_onset']] #limited to pre-stim
                        Binary_spikes1=np.zeros(math.ceil(Units_log.loc[neuron1,'stim_onset']*1000))#
                        for spike in spikes1:
                            Binary_spikes1[int(spike*1000)]=True
                        Binary_spikes1= Binary_spikes1[-1000:]# limit to only the last 1sec before stim
                        
                        
                        for j,neuron2 in enumerate(Units_log.index):
                            spikes2=Units_log.loc[neuron2,'spike_times']
                            spikes2=spikes2[spikes2<Units_log.loc[neuron2,'stim_onset']] #limited to pre-stim
                            Binary_spikes2=np.zeros(math.ceil(Units_log.loc[neuron2,'stim_onset']*1000))#
                            for spike in spikes2:
                                Binary_spikes2[int(spike*1000)]=True
                            Binary_spikes2= Binary_spikes2[-1000:]# limit to only the last 1sec before stim
                            
                           
                            corr=sp.signal.correlate(Binary_spikes1,Binary_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_spikes2))) 
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -500 to 500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            REALcorr250msecNorm= corr250msecNorm/Triangle
        
                            
                            
                            #SHUFFLE
                            #spikes2=[x[0] for x in spikes2]
                            NextTrial_spikes2=Units_log_NextTrial.loc[Units_log_NextTrial.index[j], 'spike_times']
                            NextTrial_spikes2=NextTrial_spikes2[NextTrial_spikes2<Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']] #limited to pre-stim
                            Binary_NextTrial_spikes2=np.zeros(math.ceil(Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']*1000))#
                            for spike in NextTrial_spikes2:
                                Binary_NextTrial_spikes2[int(spike*1000)]=True
                            Binary_NextTrial_spikes2= Binary_NextTrial_spikes2[-1000:]# limit to only the last 1sec before stim
                            
                            corr=sp.signal.correlate(Binary_spikes1,Binary_NextTrial_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_NextTrial_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_NextTrial_spikes2))) #ATTENTION: I'm not sure the ean for neuron2 sould include the zeros-padded?
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -500 to 500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            SHUFFLEcorr250msecNorm= corr250msecNorm/Triangle
                            
                            
                            
                            # ATTENTION EXCEPTION:#This is very important because even though the trial is being counted 
                            # as having enough time, the correlation will be zero if there are no spikes. 
                            #This messes up the normailzation by trials, especially when using small time_windows 
                            #(ie. high chance of no spikes)
                            if ((sum(Binary_spikes1)==0) | (sum(Binary_spikes2)==0) | (sum(Binary_NextTrial_spikes2)==0)): 
                                Number_of_counted_trials_PAIRMATRIX[i,j]-=1
                            else:
                                day_df.at[i,j]=day_df.loc[i,j]+REALcorr250msecNorm
                                shuffle_df.at[i,j]=shuffle_df.loc[i,j]+SHUFFLEcorr250msecNorm
    
                #    print(trial)
                
                for i in day_df.index:
                    for j in day_df.index:
                        if Number_of_counted_trials_PAIRMATRIX[i,j]>30: #do not include this session if there were fewer than 30 trials to get spikes from
                            #Normalize by number of trials
                            day_df.at[i,j]=day_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                            shuffle_df.at[i,j]=shuffle_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                        else: #Make them Zero So that you can discount them later
                            day_df.at[i,j]=np.zeros((1,1000))          
                            shuffle_df.at[i,j]=np.zeros((1,1000))
                     
                    
                    
                for i in range(Units_number):
                    for j in range(Units_number):
                        if Units_number>1:
                             if j>=i: #avoid dupictaes
                                 if ( (Names_for_plot[i][-8:-5]!=Names_for_plot[j][-8:-5]) | (Names_for_plot[i]==Names_for_plot[j])): # Dont count units on same tetrodes
                                     XCorr_df.at[Counting_pairs,'Neuron1_name']=Names_for_plot[i]
                                     XCorr_df.at[Counting_pairs,'Neuron2_name']=Names_for_plot[j]
    #                                
                                     XCorr_df.at[Counting_pairs,'Neuron1_Category']=Category_for_plot[i]
                                     XCorr_df.at[Counting_pairs,'Neuron2_Category']=Category_for_plot[j]
                                     if trial_set=='All':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_ALL1sec']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_ALL1sec']=shuffle_df.loc[i,j][0]
                                     elif trial_set=='Correct':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_CORRECT1sec']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_CORRECT1sec']=shuffle_df.loc[i,j][0]
                                     elif trial_set=='Incorrect':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_INCORRECT1sec']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_INCORRECT1sec']=shuffle_df.loc[i,j][0]
                                     elif trial_set=='Right':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHT1sec']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_RIGHT1sec']=shuffle_df.loc[i,j][0]
                                     elif trial_set=='Left':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFT1sec']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFT1sec']=shuffle_df.loc[i,j][0]
          
                                     Counting_pairs+=1

    return XCorr_df
     


###############################################################################
###############################################################################
######################### S1XCorr #################################
###############################################################################
###############################################################################

def Generate_S1XCorr_df(S1master_log)

mice=['EF0074', 'EF0076','EF0077','EF0079']  
S1XCorr_df=pd.DataFrame(columns=['Neuron1_name', 'Neuron2_name','Neuron1_LickPref','Neuron2_LickPref', 'Neuron1_Category', 'Neuron2_Category', 'Neuron1_group', 'Neuron2_group', 'XCorr-50to50_ALL1sec', 'XCorr-50to50_shuffle_ALL1sec','XCorr-50to50_RIGHT1sec', 'XCorr-50to50_shuffle_RIGHT1sec','XCorr-50to50_LEFT1sec', 'XCorr-50to50_shuffle_LEFT1sec', 'XCorr-50to50_CORRECT1sec', 'XCorr-50to50_shuffle_CORRECT1sec', 'XCorr-50to50_INCORRECT1sec', 'XCorr-50to50_shuffle_INCORRECT1sec'])
   
for trial_set in ['All','Correct', 'Incorrect', 'Right','Left']:
    Counting_pairs=0
    master_counter=0
    master_LickPref=[]
    for mouse in mice:
        print(mouse)
        mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
        day_count=0
        for day in np.unique(mouse_log['date']):   
            print(day)
            day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
            #Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num']) #these are a mess of either ints or ints trapped in arrays...

            if trial_set=='All':
                    Trial_numbers=np.unique(day_log['trial_num'])
                    #Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
            elif trial_set=='Correct':
                    Trial_numbers=np.unique(day_log[(day_log['Stim/Block/Response']=='SomCR') | (day_log['Stim/Block/Response']=='VisCR') | (day_log['Stim/Block/Response']=='SomHit') | (day_log['Stim/Block/Response']=='VisHit')]['trial_num'])
                    #Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
            elif trial_set=='Incorrect':
                    Trial_numbers=np.unique(day_log[(day_log['Stim/Block/Response']=='SomMiss') | (day_log['Stim/Block/Response']=='VisMiss') | (day_log['Stim/Block/Response']=='SomFA') | (day_log['Stim/Block/Response']=='VisFA')]['trial_num'])
                    #Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
            elif trial_set=='Right':
                    Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num'])
                    #Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)              
            elif trial_set=='Left':
                    Trial_numbers=np.unique(day_log[day_log['block_type']=='Visual']['trial_num'])
                    #Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                
            Units_number=len(np.unique(day_log['unit_name']))
                    
            Units_number=len(np.unique(day_log['unit_name']))
            
            Names_for_plot=[]
            Category_for_plot=[]
            #First loop through and get their LickPref
            for unit in np.unique(day_log['unit_name']):
                unit_log=day_log[day_log['unit_name']==unit]
                Names_for_plot.append(unit_log['unit_name'].values[0])
               
                
#                if unit_log['Category'].values[0]=='OptoTag':
#                    Category_for_plot.append('OptoTag')
#                elif unit_log['Category2'].values[0]=='OptoNetwork':
#                    Category_for_plot.append('SALT')
#                elif unit_log['Category'].values[0]=='SameTT':
#                    Category_for_plot.append('SameTT')
#                elif unit_log['Category3'].values[0]=='OptoNetwork_SameTT':
#                    Category_for_plot.append('SALT_SameTT')
#                elif unit_log['Category4'].values[0]=='InBetween':
#                    Category_for_plot.append('InBetween')
#                else:
#                    Category_for_plot.append('X')
            
            #create dat frame that will be used to store and compile spikes from all pairs and all trials
            day_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
            for row in day_df.index:
                for col in day_df.index:
                    day_df.at[row,col]=np.zeros((1,1000))
                    
            #create dat frame that will be used to store and compile spikes from all pairs and all trials
            shuffle_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
            for row in shuffle_df.index:
                for col in shuffle_df.index:
                    shuffle_df.at[row,col]=np.zeros((1,1000))
                    
            Number_of_counted_trials=   len(Trial_numbers)-1 #we will discount each trial we remove to check that we have enough and then normalize by the right amount     
            temp=np.repeat(Number_of_counted_trials,Units_number)
            Number_of_counted_trials_PAIRMATRIX=np.asarray([ temp for _ in np.arange(Units_number)])
            for m,trial in enumerate(Trial_numbers[:-1]): #here we need to enumerate to grab the next trial when doing the control (instead of shuffle) (go until the second to last because we are taking '+1')
                Units_log=day_log[day_log['trial_num']==trial] #represent one trial, on ow per neuron
                Units_log_NextTrial=day_log[day_log['trial_num']==Trial_numbers[m+1]] #will be used for shifted XCorr
                
                #Deal with error when a trial has no 'Stim_onset'               
                a=np.unique(Units_log['stim_onset'])
                if np.sum(a)==0:
                    a=0
                else:
                    a=a[0]#[0][0][0][0]
                b=np.unique( Units_log_NextTrial['stim_onset'])
                if np.sum(b)==0:
                    b=0
                else:
                    b=b[0]#[0][0][0][0]

                if ((a<1) | (b<1)):
                    Number_of_counted_trials_PAIRMATRIX[:,:]-=1
                    print(a)
                    continue
                
                
                
                for i,neuron1 in enumerate(Units_log.index):
                    spikes1=Units_log.loc[neuron1,'spike_times']
                    spikes1=spikes1[spikes1<Units_log.loc[neuron1,'stim_onset']] #limited to the 3 seconds pre-stim
                    Binary_spikes1=np.zeros(math.ceil(Units_log.loc[neuron1,'stim_onset']*1000))#
                    for spike in spikes1:
                        Binary_spikes1[int(spike*1000)]=True
                    Binary_spikes1= Binary_spikes1[-1000:]# limit to only the last 3sec before stim
                    
                    
                    for j,neuron2 in enumerate(Units_log.index):
                        spikes2=Units_log.loc[neuron2,'spike_times']
                        spikes2=spikes2[spikes2<Units_log.loc[neuron2,'stim_onset']] #limited to  pre-stim
                        Binary_spikes2=np.zeros(math.ceil(Units_log.loc[neuron2,'stim_onset']*1000))#
                        for spike in spikes2:
                            Binary_spikes2[int(spike*1000)]=True
                        Binary_spikes2= Binary_spikes2[-1000:]# limit to only the last 1sec before stim
                        
                       
                        corr=sp.signal.correlate(Binary_spikes1,Binary_spikes2)
                        corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                        
                        #Normailze by geometric mean
                        if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_spikes2)>0)):
                            corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_spikes2))) 
                        else:
                            corr250msecNorm=corr250msec
    
                        #Normalize by a triangular filter to cancel the effect of zero-padding
                        # We are using a lag from -500 to 500, each value of the 
                        # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                        # make it seconds so that the final value is /s
                        Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                        REALcorr250msecNorm= corr250msecNorm/Triangle
    
                        
                        
                        #SHUFFLE
                        #spikes2=[x[0] for x in spikes2]
                        NextTrial_spikes2=Units_log_NextTrial.loc[Units_log_NextTrial.index[j], 'spike_times']
                        NextTrial_spikes2=NextTrial_spikes2[NextTrial_spikes2<Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']] #limited to  pre-stim
                        Binary_NextTrial_spikes2=np.zeros(math.ceil(Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']*1000))#
                        for spike in NextTrial_spikes2:
                            Binary_NextTrial_spikes2[int(spike*1000)]=True
                        Binary_NextTrial_spikes2= Binary_NextTrial_spikes2[-1000:]# limit to only the last 1sec before stim
                        
                        corr=sp.signal.correlate(Binary_spikes1,Binary_NextTrial_spikes2)
                        corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                        #Normailze by geometric mean
                        if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_NextTrial_spikes2)>0)):
                            corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_NextTrial_spikes2))) 
                        else:
                            corr250msecNorm=corr250msec
    
                        #Normalize by a triangular filter to cancel the effect of zero-padding
                        # We are using a lag from -500 to 500, each value of the 
                        # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                        # make it seconds so that the final value is /s
                        Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                        SHUFFLEcorr250msecNorm= corr250msecNorm/Triangle
                        
                        
                        
                        # ATTENTION EXCEPTION:#This is very important because even though the trial is being counted 
                        # as having enough time, the correlation will be zero if there are no spikes. 
                        #This messes up the normailzation by trials, especially when using small time_windows 
                        #(ie. high chance of no spikes)
                        if ((sum(Binary_spikes1)==0) | (sum(Binary_spikes2)==0) | (sum(Binary_NextTrial_spikes2)==0)): 
                            Number_of_counted_trials_PAIRMATRIX[i,j]-=1
                        else:
                            day_df.at[i,j]=day_df.loc[i,j]+REALcorr250msecNorm
                            shuffle_df.at[i,j]=shuffle_df.loc[i,j]+SHUFFLEcorr250msecNorm

            #    print(trial)
            
            for i in day_df.index:
                for j in day_df.index:
                    if Number_of_counted_trials_PAIRMATRIX[i,j]>30: #do not include this session if there were fewer than 30 trials to get spikes from
                        #Normalize by number of trials
                        day_df.at[i,j]=day_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                        shuffle_df.at[i,j]=shuffle_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                    else: #Make them Zero So that you can discount them later
                        day_df.at[i,j]=np.zeros((1,1000))          
                        shuffle_df.at[i,j]=np.zeros((1,1000))
                 
                
                
            for i in range(Units_number):
                for j in range(Units_number):
                    if Units_number>1:
                         if j>=i: #avoid dupictaes
                             if ( (Names_for_plot[i][-8:-5]!=Names_for_plot[j][-8:-5]) | (Names_for_plot[i]==Names_for_plot[j])): # Dont count units on same tetrodes
                                 S1XCorr_df.at[Counting_pairs,'Neuron1_name']=Names_for_plot[i]
                                 S1XCorr_df.at[Counting_pairs,'Neuron2_name']=Names_for_plot[j]

#                                 S1XCorr_df.at[Counting_pairs,'Neuron1_Category']=Category_for_plot[i]
#                                 S1XCorr_df.at[Counting_pairs,'Neuron2_Category']=Category_for_plot[j]
                                 if trial_set=='All':
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_ALL1sec']=day_df.loc[i,j][0]
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_ALL1sec']=shuffle_df.loc[i,j][0]
                                 elif trial_set=='Correct':
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_CORRECT1sec']=day_df.loc[i,j][0]
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_CORRECT1sec']=shuffle_df.loc[i,j][0]
                                 elif trial_set=='Incorrect':
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_INCORRECT1sec']=day_df.loc[i,j][0]
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_INCORRECT1sec']=shuffle_df.loc[i,j][0]
                                 elif trial_set=='Right':
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHT1sec']=day_df.loc[i,j][0]
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_RIGHT1sec']=shuffle_df.loc[i,j][0]
                                 elif trial_set=='Left':
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFT1sec']=day_df.loc[i,j][0]
                                     S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFT1sec']=shuffle_df.loc[i,j][0]
                                                 
                                 Counting_pairs+=1     
    return S1XCorr_df


###############################################################################
###############################################################################
######################### Add POST to XCorr_df #################################
###############################################################################
###############################################################################
def Add_POST_to_XCorr_df(master_log, XCorr_df):
    mice=['Cl4','Cl5','Cl6','Claustrum18','Claustrum23','Claustrum25','Claustrum31','Claustrum32','Claustrum37']
    #Now XCorr   
    #XCorr_df=pd.DataFrame(columns=['Neuron1_name', 'Neuron2_name','Neuron1_LickPref','Neuron2_LickPref', 'Neuron1_Category', 'Neuron2_Category', 'Neuron1_group', 'Neuron2_group', 'XCorr-50to50_ALL1sec', 'XCorr-50to50_shuffle_ALL1sec','XCorr-50to50_RIGHT1sec', 'XCorr-50to50_shuffle_RIGHT1sec','XCorr-50to50_LEFT1sec', 'XCorr-50to50_shuffle_LEFT1sec', 'XCorr-50to50_CORRECT1sec', 'XCorr-50to50_shuffle_CORRECT1sec', 'XCorr-50to50_INCORRECT1sec', 'XCorr-50to50_shuffle_INCORRECT1sec'])
    
    ##                                 
    ##                                 
    Name_df=['XCorr-50to50_ALL1secPOST','XCorr-50to50_shuffle_ALL1secPOST',
             'XCorr-50to50_RIGHTHIT1secPOST','XCorr-50to50_shuffle_RIGHTHIT1secPOST',
             'XCorr-50to50_LEFTHIT1secPOST','XCorr-50to50_shuffle_LEFTHIT1secPOST', 
             'XCorr-50to50_RIGHTCR1secPOST','XCorr-50to50_shuffle_RIGHTCR1secPOST',
             'XCorr-50to50_LEFTCR1secPOST','XCorr-50to50_shuffle_LEFTCR1secPOST',
             'Geometric mean FR ALLPOST',
             'Geometric mean FR RIGHTHITPOST', 'Geometric mean FR LEFTHITPOST',
             'Geometric mean FR RIGHTCRPOST', 'Geometric mean FR LEFTCRPOST']
    for i,name in enumerate(Name_df):   
        XCorr_df[Name_df[i]]=np.zeros(len(XCorr_df))
        XCorr_df[Name_df[i]]=XCorr_df[Name_df[i]].astype(object)
    #end of preprocessing                                 
    ##                                 
    
    for trial_set in ['All','RightHit','LeftHit','RightCR', 'LeftCR']:
        Counting_pairs=0
        master_counter=0
        master_LickPref=[]
        Geometric_mean_FR=[]
        for mouse in mice:
            print(mouse)
            mouse_log=master_log.loc[np.equal(master_log['mouse_name'], mouse)]
            day_count=0
            for day in np.unique(mouse_log['date']):   
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                #Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num']) #these are a mess of either ints or ints trapped in arrays...
    
                if trial_set=='All':
                        Trial_numbers=np.unique(day_log['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='RightHit':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomHit']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisHit']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='LeftHit':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisHit']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomHit']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)        
                elif trial_set=='RightCR':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomCR']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisCR']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='LeftCR':
                    if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisCR']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                    elif mouse in REVERSED:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomCR']['trial_num'])
                        Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                Units_number=len(np.unique(day_log['unit_name']))
                        
                Units_number=len(np.unique(day_log['unit_name']))
                
                Lick_prefs=[]
                Names_for_plot=[]
                Category_for_plot=[]
                #First loop through and get their LickPref
                for unit in np.unique(day_log['unit_name']):
                    unit_log=day_log[day_log['unit_name']==unit]
                    Names_for_plot.append(unit_log['unit_name'].values[0])
    #                Lick_prefs.append(unit_log['LickPref'].values[0])
    
                    if unit_log['Category'].values[0]=='OptoTag':
                        Category_for_plot.append('OptoTag')
                    elif unit_log['Category2'].values[0]=='OptoNetwork':
                        Category_for_plot.append('SALT')
                    elif unit_log['Category'].values[0]=='SameTT':
                        Category_for_plot.append('SameTT')
                    elif unit_log['Category3'].values[0]=='OptoNetwork_SameTT':
                        Category_for_plot.append('SALT_SameTT')
                    elif unit_log['Category4'].values[0]=='InBetween':
                        Category_for_plot.append('InBetween')
                    else:
                        Category_for_plot.append('X')
                
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                day_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in day_df.index:
                    for col in day_df.index:
                        day_df.at[row,col]=np.zeros((1,1000))
                        
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                shuffle_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in shuffle_df.index:
                    for col in shuffle_df.index:
                        shuffle_df.at[row,col]=np.zeros((1,1000))
                
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                FR_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in FR_df.index:
                    for col in FR_df.index:
                        FR_df.at[row,col]=np.zeros((1,2))
                        
                Number_of_counted_trials=   len(Trial_numbers)-1 #we will discount each trial we remove to check that we have enough and then normalize by the right amount     
                temp=np.repeat(Number_of_counted_trials,Units_number)
                Number_of_counted_trials_PAIRMATRIX=np.asarray([ temp for _ in np.arange(Units_number)])
                for m,trial in enumerate(Trial_numbers[:-1]): #here we need to enumerate to grab the next trial when doing the control (instead of shuffle) (go until the second to last because we are taking '+1')
                    Units_log=day_log[day_log['trial_num']==trial] #represent one trial, on ow per neuron
                    Units_log_NextTrial=day_log[day_log['trial_num']==Trial_numbers[m+1]] #will be used for shifted XCorr
                    
                    #Deal with error when a trial has no 'Stim_onset'               
                    a=np.unique(Units_log['stim_onset'])
                    if np.sum(a)==0:
                        a=0
                    else:
                        a=a[0][0][0][0][0]
                    b=np.unique( Units_log_NextTrial['stim_onset'])
                    if np.sum(b)==0:
                        b=0
                    else:
                        b=b[0][0][0][0][0]
    
                    if ((a<1) | (b<1)):
                        Number_of_counted_trials_PAIRMATRIX[:,:]-=1
                        print(a)
                        continue
                    
                    
                    
                    for i,neuron1 in enumerate(Units_log.index):
                        spikes1=Units_log.loc[neuron1,'spike_times']
                        spikes1=spikes1[spikes1<(Units_log.loc[neuron1,'stim_onset']+1)] #get spikes until Stim+1sec
                        Binary_spikes1=np.zeros(math.ceil((Units_log.loc[neuron1,'stim_onset']+1)*1000))#
                        for spike in spikes1:
                            Binary_spikes1[int(spike*1000)]=True
                        Binary_spikes1= Binary_spikes1[-1000:]# get the last seconds, which is the 1 sec post stim
                        SpikeRate1=sum(Binary_spikes1)
                        
                        
                        for j,neuron2 in enumerate(Units_log.index):
                            spikes2=Units_log.loc[neuron2,'spike_times']
                            spikes2=spikes2[spikes2<(Units_log.loc[neuron2,'stim_onset']+1)] #get spikes until Stim+1sec
                            Binary_spikes2=np.zeros(math.ceil((Units_log.loc[neuron2,'stim_onset']+1)*1000))#
                            for spike in spikes2:
                                Binary_spikes2[int(spike*1000)]=True
                            Binary_spikes2= Binary_spikes2[-1000:]# get the last seconds, which is the 1 sec post stim
                            SpikeRate2=sum(Binary_spikes2)
                           
                            corr=sp.signal.correlate(Binary_spikes1,Binary_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_spikes2))) 
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -1500 to 1500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            REALcorr250msecNorm= corr250msecNorm/Triangle
        
                            
                            
                            #SHUFFLE
                            #spikes2=[x[0] for x in spikes2]
                            NextTrial_spikes2=Units_log_NextTrial.loc[Units_log_NextTrial.index[j], 'spike_times']
                            NextTrial_spikes2=NextTrial_spikes2[NextTrial_spikes2<(Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']+1)] #get spikes until Stim+1sec
                            Binary_NextTrial_spikes2=np.zeros(math.ceil((Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']+1)*1000))#
                            for spike in NextTrial_spikes2:
                                Binary_NextTrial_spikes2[int(spike*1000)]=True
                            Binary_NextTrial_spikes2= Binary_NextTrial_spikes2[-1000:]# get the last seconds, which is the 1 sec post stim
                            
                            corr=sp.signal.correlate(Binary_spikes1,Binary_NextTrial_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_NextTrial_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_NextTrial_spikes2))) #ATTENTION: I'm not sure the ean for neuron2 sould include the zeros-padded?
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -1500 to 1500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            SHUFFLEcorr250msecNorm= corr250msecNorm/Triangle
                            
                            
                            
                            # ATTENTION EXCEPTION:#This is very important because even though the trial is being counted 
                            # as having enough time, the correlation will be zero if there are no spikes. 
                            #This messes up the normailzation by trials, especially when using small time_windows 
                            #(ie. high chance of no spikes)
                            if ((sum(Binary_spikes1)==0) | (sum(Binary_spikes2)==0) | (sum(Binary_NextTrial_spikes2)==0)): 
                                Number_of_counted_trials_PAIRMATRIX[i,j]-=1
                            else:
                                day_df.at[i,j]=day_df.loc[i,j]+REALcorr250msecNorm
                                shuffle_df.at[i,j]=shuffle_df.loc[i,j]+SHUFFLEcorr250msecNorm
                                FR_df.at[i,j]=FR_df.loc[i,j]+np.array([SpikeRate1, SpikeRate2])
    
                #    print(trial)
                
                for i in day_df.index:
                    for j in day_df.index:
                        if Number_of_counted_trials_PAIRMATRIX[i,j]>30: #do not include this session if there were fewer than 30 trials to get spikes from
                            #Normalize by number of trials
                            day_df.at[i,j]=day_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                            shuffle_df.at[i,j]=shuffle_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                        else: #Make them Zero So that you can discount them later
                            day_df.at[i,j]=np.zeros((1,1000))          
                            shuffle_df.at[i,j]=np.zeros((1,1000))
                            FR_df.at[i,j][0][0]=0 
                            FR_df.at[i,j][0][1]=0 
                     
                    
                    
                for i in range(Units_number):
                    for j in range(Units_number):
                        if Units_number>1:
                             if j>=i: #avoid dupictaes
                                 if ( (Names_for_plot[i][-8:-5]!=Names_for_plot[j][-8:-5]) | (Names_for_plot[i]==Names_for_plot[j])): # Dont count units on same tetrodes
                                     XCorr_df.at[Counting_pairs,'Neuron1_name']=Names_for_plot[i]
                                     XCorr_df.at[Counting_pairs,'Neuron2_name']=Names_for_plot[j]

                                     XCorr_df.at[Counting_pairs,'Neuron1_Category']=Category_for_plot[i]
                                     XCorr_df.at[Counting_pairs,'Neuron2_Category']=Category_for_plot[j]
                                     if trial_set=='All':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_ALL1secPOST']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_ALL1secPOST']=shuffle_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,'Geometric mean FR ALLPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
                                     
                                     elif trial_set=='RightHit':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHTHIT1secPOST']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,  'XCorr-50to50_shuffle_RIGHTHIT1secPOST']=shuffle_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,'Geometric mean FR RIGHTHITPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='LeftHit':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFTHIT1secPOST']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFTHIT1secPOST']=shuffle_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,'Geometric mean FR LEFTHITPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='RightCR':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHTCR1secPOST']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_RIGHTCR1secPOST']=shuffle_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,'Geometric mean FR RIGHTCRPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='LeftCR':
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFTCR1secPOST']=day_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFTCR1secPOST']=shuffle_df.loc[i,j][0]
                                         XCorr_df.at[Counting_pairs,'Geometric mean FR LEFTCRPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]

                                     Counting_pairs+=1
    return XCorr_df




###############################################################################
###############################################################################
######################### Add POST to S1XCorr_df #################################
###############################################################################
###############################################################################
def Add_POST_to_S1XCorr_df(master_log, XCorr_df):
    mice=['EF0074', 'EF0076','EF0077','EF0079']  
       
    #XCorr_df=pd.DataFrame(columns=['Neuron1_name', 'Neuron2_name','Neuron1_LickPref','Neuron2_LickPref', 'Neuron1_Category', 'Neuron2_Category', 'Neuron1_group', 'Neuron2_group', 'XCorr-50to50_ALL1sec', 'XCorr-50to50_shuffle_ALL1sec','XCorr-50to50_RIGHT1sec', 'XCorr-50to50_shuffle_RIGHT1sec','XCorr-50to50_LEFT1sec', 'XCorr-50to50_shuffle_LEFT1sec', 'XCorr-50to50_CORRECT1sec', 'XCorr-50to50_shuffle_CORRECT1sec', 'XCorr-50to50_INCORRECT1sec', 'XCorr-50to50_shuffle_INCORRECT1sec'])
    
    ##                                 
    ##                                 
    Name_df=['XCorr-50to50_ALL1secPOST','XCorr-50to50_shuffle_ALL1secPOST',
             'XCorr-50to50_RIGHTHIT1secPOST','XCorr-50to50_shuffle_RIGHTHIT1secPOST',
             'XCorr-50to50_LEFTHIT1secPOST','XCorr-50to50_shuffle_LEFTHIT1secPOST', 
             'XCorr-50to50_RIGHTCR1secPOST','XCorr-50to50_shuffle_RIGHTCR1secPOST',
             'XCorr-50to50_LEFTCR1secPOST','XCorr-50to50_shuffle_LEFTCR1secPOST',
             'Geometric mean FR ALLPOST',
             'Geometric mean FR RIGHTHITPOST', 'Geometric mean FR LEFTHITPOST',
             'Geometric mean FR RIGHTCRPOST', 'Geometric mean FR LEFTCRPOST']
    for i,name in enumerate(Name_df):   
        S1XCorr_df[Name_df[i]]=np.zeros(len(S1XCorr_df))
        S1XCorr_df[Name_df[i]]=S1XCorr_df[Name_df[i]].astype(object)
    #end of preprocessing                                 
    ##                                 
    
    for trial_set in ['All','RightHit','LeftHit','RightCR', 'LeftCR']:
        Counting_pairs=0
        master_counter=0
        master_LickPref=[]
        Geometric_mean_FR=[]
        for mouse in mice:
            print(mouse)
            mouse_log=S1master_log.loc[np.equal(S1master_log['mouse_name'], mouse)]
            day_count=0
            for day in np.unique(mouse_log['date']):   
                print(day)
                day_log=mouse_log.loc[np.equal(mouse_log['date'], day)]
                #Trial_numbers=np.unique(day_log[day_log['block_type']=='Whisker']['trial_num']) #these are a mess of either ints or ints trapped in arrays...
    
                if trial_set=='All':
                        Trial_numbers=np.unique(day_log['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='RightHit':
    #                if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomHit']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
    #                elif mouse in REVERSED:
    #                    Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisHit']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='LeftHit':
    #                if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisHit']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
    #                elif mouse in REVERSED:
    #                    Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomHit']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)        
                elif trial_set=='RightCR':
    #                if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomCR']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
    #                elif mouse in REVERSED:
    #                    Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisCR']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                elif trial_set=='LeftCR':
    #                if mouse in ORIGINAL:
                        Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='VisCR']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
    #                elif mouse in REVERSED:
    #                    Trial_numbers=np.unique(day_log[day_log['Stim/Block/Response']=='SomCR']['trial_num'])
    #                    Trial_numbers=[x if isinstance(x,int) else x[0][0][0][0] for x in Trial_numbers]#...cleaned up :)
                Units_number=len(np.unique(day_log['unit_name']))
                        
                
                Lick_prefs=[]
                Names_for_plot=[]
                Category_for_plot=[]
                #First loop through and get their LickPref
                for unit in np.unique(day_log['unit_name']):
                    unit_log=day_log[day_log['unit_name']==unit]
                    Names_for_plot.append(unit_log['unit_name'].values[0])
    #                Lick_prefs.append(unit_log['LickPref'].values[0])
    
    #                if unit_log['Category'].values[0]=='OptoTag':
    #                    Category_for_plot.append('OptoTag')
    #                elif unit_log['Category2'].values[0]=='OptoNetwork':
    #                    Category_for_plot.append('SALT')
    #                elif unit_log['Category'].values[0]=='SameTT':
    #                    Category_for_plot.append('SameTT')
    #                elif unit_log['Category3'].values[0]=='OptoNetwork_SameTT':
    #                    Category_for_plot.append('SALT_SameTT')
    #                elif unit_log['Category4'].values[0]=='InBetween':
    #                    Category_for_plot.append('InBetween')
    #                else:
    #                    Category_for_plot.append('X')
                
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                day_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in day_df.index:
                    for col in day_df.index:
                        day_df.at[row,col]=np.zeros((1,1000))
                        
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                shuffle_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in shuffle_df.index:
                    for col in shuffle_df.index:
                        shuffle_df.at[row,col]=np.zeros((1,1000))
                
                #create dat frame that will be used to store and compile spikes from all pairs and all trials
                FR_df=pd.DataFrame(index=np.arange(Units_number), columns=np.arange(Units_number))
                for row in FR_df.index:
                    for col in FR_df.index:
                        FR_df.at[row,col]=np.zeros((1,2))
                        
                Number_of_counted_trials=   len(Trial_numbers)-1 #we will discount each trial we remove to check that we have enough and then normalize by the right amount     
                temp=np.repeat(Number_of_counted_trials,Units_number)
                Number_of_counted_trials_PAIRMATRIX=np.asarray([ temp for _ in np.arange(Units_number)])
                for m,trial in enumerate(Trial_numbers[:-1]): #here we need to enumerate to grab the next trial when doing the control (instead of shuffle) (go until the second to last because we are taking '+1')
                    Units_log=day_log[day_log['trial_num']==trial] #represent one trial, on ow per neuron
                    Units_log_NextTrial=day_log[day_log['trial_num']==Trial_numbers[m+1]] #will be used for shifted XCorr
                    
                    #Deal with error when a trial has no 'Stim_onset'               
                    a=np.unique(Units_log['stim_onset'])
                    if np.sum(a)==0:
                        a=0
                    else:
                        a=a[0]
                    b=np.unique( Units_log_NextTrial['stim_onset'])
                    if np.sum(b)==0:
                        b=0
                    else:
                        b=b[0]
    
                    if ((a<1) | (b<1)):
                        Number_of_counted_trials_PAIRMATRIX[:,:]-=1
                        print(a)
                        continue
                    
                    
                    
                    for i,neuron1 in enumerate(Units_log.index):
                        spikes1=Units_log.loc[neuron1,'spike_times']
                        spikes1=spikes1[spikes1<(Units_log.loc[neuron1,'stim_onset']+1)] #get spikes until Stim+1sec
                        Binary_spikes1=np.zeros(math.ceil((Units_log.loc[neuron1,'stim_onset']+1)*1000))#
                        for spike in spikes1:
                            Binary_spikes1[int(spike*1000)]=True
                        Binary_spikes1= Binary_spikes1[-1000:]# get the last seconds, which is the 1 sec post stim
                        SpikeRate1=sum(Binary_spikes1)
                        
                        
                        for j,neuron2 in enumerate(Units_log.index):
                            spikes2=Units_log.loc[neuron2,'spike_times']
                            spikes2=spikes2[spikes2<(Units_log.loc[neuron2,'stim_onset']+1)] #get spikes until Stim+1sec
                            Binary_spikes2=np.zeros(math.ceil((Units_log.loc[neuron2,'stim_onset']+1)*1000))#
                            for spike in spikes2:
                                Binary_spikes2[int(spike*1000)]=True
                            Binary_spikes2= Binary_spikes2[-1000:]# get the last seconds, which is the 1 sec post stim
                            SpikeRate2=sum(Binary_spikes2)
                           
                            corr=sp.signal.correlate(Binary_spikes1,Binary_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_spikes2))) 
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -500 to 500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            REALcorr250msecNorm= corr250msecNorm/Triangle
        
                            
                            
                            #SHUFFLE
                            #spikes2=[x[0] for x in spikes2]
                            NextTrial_spikes2=Units_log_NextTrial.loc[Units_log_NextTrial.index[j], 'spike_times']
                            NextTrial_spikes2=NextTrial_spikes2[NextTrial_spikes2<(Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']+1)] #get spikes until Stim+1sec
                            Binary_NextTrial_spikes2=np.zeros(math.ceil((Units_log_NextTrial.loc[Units_log_NextTrial.index[j],'stim_onset']+1)*1000))#
                            for spike in NextTrial_spikes2:
                                Binary_NextTrial_spikes2[int(spike*1000)]=True
                            Binary_NextTrial_spikes2= Binary_NextTrial_spikes2[-1000:]# get the last seconds, which is the 1 sec post stim
                            
                            corr=sp.signal.correlate(Binary_spikes1,Binary_NextTrial_spikes2)
                            corr250msec=corr[len( Binary_spikes1)-500: len(Binary_spikes1)+500]
                            #Normailze by geometric mean
                            if ((np.mean(Binary_spikes1)>0) & (np.mean(Binary_NextTrial_spikes2)>0)):
                                corr250msecNorm=corr250msec/(np.sqrt(np.sum(Binary_spikes1) * np.sum(Binary_NextTrial_spikes2))) #ATTENTION: I'm not sure the ean for neuron2 sould include the zeros-padded?
                            else:
                                corr250msecNorm=corr250msec
        
                            #Normalize by a triangular filter to cancel the effect of zero-padding
                            # We are using a lag from -500 to 500, each value of the 
                            # correlogram at lag X must be normalized by (TotalTime - abs(lag))
                            # make it seconds so that the final value is /s
                            Triangle=[1-(abs(lag)/1000) for lag in np.arange(-500,500)]
                            SHUFFLEcorr250msecNorm= corr250msecNorm/Triangle
                            
                            
                            
                            # ATTENTION EXCEPTION:#This is very important because even though the trial is being counted 
                            # as having enough time, the correlation will be zero if there are no spikes. 
                            #This messes up the normailzation by trials, especially when using small time_windows 
                            #(ie. high chance of no spikes)
                            if ((sum(Binary_spikes1)==0) | (sum(Binary_spikes2)==0) | (sum(Binary_NextTrial_spikes2)==0)): 
                                Number_of_counted_trials_PAIRMATRIX[i,j]-=1
                            else:
                                day_df.at[i,j]=day_df.loc[i,j]+REALcorr250msecNorm
                                shuffle_df.at[i,j]=shuffle_df.loc[i,j]+SHUFFLEcorr250msecNorm
                                FR_df.at[i,j]=FR_df.loc[i,j]+np.array([SpikeRate1, SpikeRate2])
    
                #    print(trial)
                
                for i in day_df.index:
                    for j in day_df.index:
                        if Number_of_counted_trials_PAIRMATRIX[i,j]>30: #do not include this session if there were fewer than 30 trials to get spikes from
                            #Normalize by number of trials
                            day_df.at[i,j]=day_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                            shuffle_df.at[i,j]=shuffle_df.loc[i,j]/Number_of_counted_trials_PAIRMATRIX[i,j]
                        else: #Make them Zero So that you can discount them later
                            day_df.at[i,j]=np.zeros((1,1000))          
                            shuffle_df.at[i,j]=np.zeros((1,1000))
                            FR_df.at[i,j][0][0]=0 
                            FR_df.at[i,j][0][1]=0 
                     
                    
                    
                for i in range(Units_number):
                    for j in range(Units_number):
                        if Units_number>1:
                             if j>=i: #avoid dupictaes
                                 if ( (Names_for_plot[i][-8:-5]!=Names_for_plot[j][-8:-5]) | (Names_for_plot[i]==Names_for_plot[j])): # Dont count units on same tetrodes
                                     S1XCorr_df.at[Counting_pairs,'Neuron1_name']=Names_for_plot[i]
                                     S1XCorr_df.at[Counting_pairs,'Neuron2_name']=Names_for_plot[j]

    #                                 S1XCorr_df.at[Counting_pairs,'Neuron1_Category']=Category_for_plot[i]
    #                                 S1XCorr_df.at[Counting_pairs,'Neuron2_Category']=Category_for_plot[j]
                                     if trial_set=='All':
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_ALL1secPOST']=day_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_ALL1secPOST']=shuffle_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,'Geometric mean FR ALLPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
                                     
                                     elif trial_set=='RightHit':
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHTHIT1secPOST']=day_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,  'XCorr-50to50_shuffle_RIGHTHIT1secPOST']=shuffle_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,'Geometric mean FR RIGHTHITPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='LeftHit':
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFTHIT1secPOST']=day_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFTHIT1secPOST']=shuffle_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,'Geometric mean FR LEFTHITPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='RightCR':
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_RIGHTCR1secPOST']=day_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_RIGHTCR1secPOST']=shuffle_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,'Geometric mean FR RIGHTCRPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]
    
                                     elif trial_set=='LeftCR':
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_LEFTCR1secPOST']=day_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs, 'XCorr-50to50_shuffle_LEFTCR1secPOST']=shuffle_df.loc[i,j][0]
                                         S1XCorr_df.at[Counting_pairs,'Geometric mean FR LEFTCRPOST']=np.sqrt(FR_df.loc[i,j][0][0] * FR_df.loc[i,j][0][1])/ Number_of_counted_trials_PAIRMATRIX[i,j]

                                     Counting_pairs+=1
    return S1XCorr_df