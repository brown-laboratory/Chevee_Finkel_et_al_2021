def Get_formatted_OptoLog():
    # This function is a copy of 'EF_OptoAnalysis_BrownLabRig.py' from which I removed 
    # the plotting and some extr analisys. Here we only use it to get the 'all_opto_metrics_df' 
    # dataframe that has the evoked spikes and it is then used in another script (which calls 
    # this function)) to calculate and plot the first spike jitter
    
    
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
    ##.............................................................................
    ##.............................................................................
    #task_data_files = glob.glob("Log_*")
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
    
    for file_num in tnrange(len(glob.glob("Opto_log_*"))):
         #mat = sp.io.loadmat(task_data_files[file_num])
         mat2 = sp.io.loadmat(opto_data_files[file_num])
        
         #log = mat['log']
         log2 = mat2['optoTable']
    
         #indv_log_df = pd.DataFrame(log, columns = column_names1)
    
         #log_df = pd.concat([log_df,indv_log_df])
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
    
    import matplotlib.patches as patches
    from matplotlib import gridspec
    # font = {'family' : 'sans-serif',
    #         'weight' : 'normal',
    #         'size'   : 18}
    
    # mpl.rc('font', **font)
    # mpl.rc('xtick', labelsize=16) 
    # mpl.rc('ytick', labelsize=16)
    # mpl.rc('axes', labelsize=18)
    
    
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
    
    #     rs1 = fig.add_subplot(gs[2])
    #     rs2 = fig.add_subplot(gs[3])
    #     rs3 = fig.add_subplot(gs[8])
    #     rs4 = fig.add_subplot(gs[9])
        
    #     ax1 = fig.add_subplot(gs[0], sharex = rs1)
    #     ax2 = fig.add_subplot(gs[1], sharex = rs2)
    #     ax3 = fig.add_subplot(gs[6], sharex = rs3)
    #     ax4 = fig.add_subplot(gs[7], sharex = rs4)
    
    #     rasters = [rs1, rs2, rs3, rs4]
    #     opto_axes = [ax1, ax2, ax3, ax4]
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
    ##............................................................................
    
    ##............................................................................
    
    def find_latencies(frequency_row):
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
    ##.............................................................................
    ##.............................................................................
    #opto_spikes_df['mouse_name'] = opto_spikes_df['mouse_name'].apply(lambda x: x[0:2]+x[-1])
    #opto_waves_df['mouse_name'] = opto_waves_df['mouse_name'].apply(lambda x: x[0:2]+x[-1])
    
    
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
        latency_metrics = opto_metrics_df.apply(lambda y: find_latencies(y), axis = 1) #find_latencies function defined above
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
    return  all_opto_metrics_df
