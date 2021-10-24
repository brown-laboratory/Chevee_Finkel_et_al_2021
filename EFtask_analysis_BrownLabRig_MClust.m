
clear
MouseName='Claustrum39';
Date='20191216';
directory=cd;
cat_session5MC(directory);
MClust
%% Read intan file and store variables
OptoStimStarts=[];
OptoStimEnds=[];
RightRewardStarts=[];
RightRewardEnds=[];
LeftRewardStarts=[];
LeftRewardEnds=[];
RightLickStarts=[];
RightLickEnds=[];
LeftLickStarts=[];
LeftLickEnds=[];
WhiskerStimStarts=[];
VisualStimStarts=[];
VisualStimEnds=[];
WhiskerBlockStarts=[];
VisualBlockStarts=[];
TrialStartTimes=[];
Full_board_dig_in_data=[];

fnall=arrayfun(@(x) x.name(1:(end)), dir('*.rhd'), 'UniformOutput', false);

for FileNum= 1 :numel(fnall)%for loop below extracts and concatenates digital inputs (could be generated when running kilosort to gain time)
    file= cell2str(fnall(FileNum)); % because of errors in cell2str, you have to remove ' and '; when you want to call the string containing just the file name
    if FileNum==1
        [notes,frequency_parameters,amplifier_channels,amplifier_data,t_amplifier,...
            aux_input_channels,aux_input_data,t_aux_input, t_dig, board_dig_in_data, board_adc_data]= Intan.read_Intan_RHD2000_file_noUI(file(3:end-3));
        numChan=numel(amplifier_channels);
        fs=frequency_parameters.amplifier_sample_rate;
    else
        [notes,frequency_parameters,amplifier_channels,amplifier_data,t_amplifier,...
            aux_input_channels,aux_input_data,t_aux_input, t_dig, board_dig_in_data, board_adc_data]= Intan.read_Intan_RHD2000_file_noUI(file(3:end-3));
    end
    
    %Create the full array because I use i to get block info below
    Full_board_dig_in_data=[ Full_board_dig_in_data, board_dig_in_data];
    
    %Need this quick if clause to define what Block we start in
    if FileNum==1
        if board_dig_in_data(9,1)==0
            WhiskerBlockStarts=t_dig(1);
        elseif board_dig_in_data(9,1)==1
            VisualBlockStarts=t_dig(1);
        end
    end
    
    
    for i=1:length(t_dig)-1 %rea each digital input and build the vectors
        
        %OptoPulses: digital input 8
        if board_dig_in_data(8,i)<board_dig_in_data(8,i+1)
            OptoStimStarts=[OptoStimStarts, t_dig(i)];
        elseif  board_dig_in_data(8,i)>board_dig_in_data(8,i+1)
            OptoStimEnds=[OptoStimEnds, t_dig(i)];
        end
        
        %Right reward detect: digital input 3
        if board_dig_in_data(3,i)<board_dig_in_data(3,i+1)
            RightRewardStarts=[RightRewardStarts, t_dig(i)];
        elseif  board_dig_in_data(3,i)>board_dig_in_data(3,i+1)
            RightRewardEnds=[RightRewardEnds, t_dig(i)];
        end
        
        %Left reward detect: digital input 4
        if board_dig_in_data(4,i)<board_dig_in_data(4,i+1)
            LeftRewardStarts=[LeftRewardStarts, t_dig(i)];
        elseif  board_dig_in_data(4,i)>board_dig_in_data(4,i+1)
            LeftRewardEnds=[ LeftRewardEnds, t_dig(i)];
        end
        
        %make Right lick Detect: digital input 1
        if board_dig_in_data(1,i)<board_dig_in_data(1,i+1)
            RightLickStarts=[RightLickStarts, t_dig(i)];
        elseif  board_dig_in_data(1,i)>board_dig_in_data(1,i+1)
            RightLickEnds=[RightLickEnds, t_dig(i)];
        end
        
        %Left lick detect: digital input 2
        if board_dig_in_data(2,i)<board_dig_in_data(2,i+1)
            LeftLickStarts=[LeftLickStarts, t_dig(i)];
        elseif  board_dig_in_data(2,i)>board_dig_in_data(2,i+1)
            LeftLickEnds=[LeftLickEnds, t_dig(i)];
        end
        
        % Whisker stim ON: digital input 5 (to get the OFF, you need to extract/write down the delay from the arduino that controls opto because the signal transits throuhg there)
        if board_dig_in_data(5,i)<board_dig_in_data(5,i+1)
            WhiskerStimStarts=[ WhiskerStimStarts, t_dig(i)];
            
        end
        
        %Visual stim: digital input 6
        if board_dig_in_data(6,i)<board_dig_in_data(6,i+1)
            VisualStimStarts=[VisualStimStarts, t_dig(i)];
        elseif  board_dig_in_data(6,i)>board_dig_in_data(6,i+1)
            VisualStimEnds=[VisualStimEnds, t_dig(i)];
        end
        
        %Block: digital input 9 (1 means visual, 0 means whisker)
        
        if board_dig_in_data(9,i)<board_dig_in_data(9,i+1)
            VisualBlockStarts=[VisualBlockStarts, t_dig(i)];
        elseif  board_dig_in_data(9,i)>board_dig_in_data(9,i+1)
            WhiskerBlockStarts=[WhiskerBlockStarts, t_dig(i)];
        end
        
        %Trial Start: digital input 1 (it is always HIGH, drops to LOW at the transition between trials)
        
        if board_dig_in_data(7,i)<board_dig_in_data(7,i+1)
            TrialStartTimes=[TrialStartTimes, t_dig(i)];
        end
    end
    
    FileNum
    
end

%% Build Log, one row per trial
%Column 1: Mouse name - 'Cl5; string
%Column 2: Date - '06-05-17' string
%Column 3: Block - 'Visual' or 'Whisker' string
%Column 4: Stim - 'Stim_Som_NoCue' or 'Stim_Vis_NoCue' string
%Column 5: Type of whisker stim - XX (this info is not extractable from the set up in the brown lab, need to record separately)
%Column 6: Type of visual stim - XX (this info is not extractable from the set up in brown lab, need to record separately
%Column 7: Response - noLick:0 lickRight:1 lickLeft:2 double
%Column 8: Trial number 1x1 cell
%Column 9: Stim onset 1x1 cell
%Column 10: Stim offset 1x1 cell (this metric would only be obtainable for visual, the length of the whiker is determined by the opto-arduino and not cpmmunicated to intan)
%Column 11: Lick Right 1x1 cell - {[1.23, 5.344]}, if empty {[]}
%Column 12: Lick Left 1x1 cell - {[1.23, 5.344]}, if empty {[]}
%Column 13: Spike times double [0.0507000000000000;0.393600000000000;3.70970000000000;4.00960000000000;4.26190000000000;4.89970000000000;5.57830000000000;6.50420000000000;6.68280000000000;8.62160000000000]
%Column 14: Cluster name string 'TT8clst1'
%
%Row organization:
%trial 1 ... ... ... TT1clst1
%trial 1 ... ... ... TT1clst2
%...
%trial 300 ... ... ... TT8clst6

%Log=cell(length(TrialStartTimes*[[[Number of units]]]), 14); % this would create an empty cell with the right size
Log=cell(length(TrialStartTimes), 14);

%these are the temporary attributes of the current block, we will acquire
%them and build the row at the end

trialNum=0;
%for i=1:length(board_dig_in-data); % we will loop through each 1/30000th of a second and update the current conditions as we go
for i=1:length(TrialStartTimes)-1
    trialNum=i;
    Log{trialNum,8}={trialNum};
    
    Log{trialNum,1}= MouseName;
    Log{trialNum, 2}= Date;
    
    %get the boudaries of the trial
    BoundariesSec=[TrialStartTimes(trialNum), TrialStartTimes(trialNum+1)]; %start and stop of trial in seconds
    T1=TrialStartTimes(trialNum)*30000; %index of start in 30000Hz timestamps
    T2=TrialStartTimes(trialNum+1)*30000; %index of stop in 30000Hz timestamps
    BoundariesTS=[ceil(T1):ceil(T2)]; %Start and Stop of trial in TimeStamp indices
    
    %Block
    if sum(Full_board_dig_in_data(9,BoundariesTS))==0 % this input is LOW during whisker blocks and HIGH during visual blocks
        Block='Whisker';
    else
        Block='Visual';
    end
    Log{trialNum,3}=Block;
    
    %Stim and StimOnset
    if ~isempty(VisualStimStarts(VisualStimStarts>BoundariesSec(1) & VisualStimStarts<BoundariesSec(2))) % if there is a visual stim within the trial
        Stim='Stim_Vis_NoCue';
        StimOnset=VisualStimStarts(VisualStimStarts>BoundariesSec(1) & VisualStimStarts<BoundariesSec(2))-TrialStartTimes(trialNum); %get the time of that stim and align to start of trial
    elseif ~isempty(WhiskerStimStarts(WhiskerStimStarts>BoundariesSec(1) & WhiskerStimStarts<BoundariesSec(2))) % if there is a whisker stim within the trial
        Stim='Stim_Som_NoCue';
        StimOnset=WhiskerStimStarts(WhiskerStimStarts>BoundariesSec(1) & WhiskerStimStarts<BoundariesSec(2))-TrialStartTimes(trialNum); %get the time of that stim and align to start of trial
    else
        Stim='NAN';
        StimOnset=[];
    end
    
    if length(StimOnset)>1 % this clause is to deal with the current [roblematic scenarios when the vis stim doesn't turn off....
        StimOnset=StimOnset(1);
    end
    
    Log{trialNum,4}=Stim;
    Log{trialNum,9}={StimOnset};
    
    %RightLick
    RightLick=RightLickStarts(RightLickStarts>BoundariesSec(1) & RightLickStarts<BoundariesSec(1)+60); %get the right licks that are within 1min of the trial start (that way you definitely get the licks after the stim has been delivered but you don't miss them if trial was terminated by manual reward delivery, for example)
    if isempty(RightLick) %if no Licks, make it 0 to avoid problems of not being to index an empty array in if statements below
        RightLick=0;
    else
        RightLick=RightLick-TrialStartTimes(trialNum); %convert to aligned to trial onset
        RightLick=RightLick'; %transpose (Eric's is a column)
    end
    Log{trialNum,11}={RightLick};
    %LeftLick
    LeftLick=LeftLickStarts(LeftLickStarts>BoundariesSec(1) & LeftLickStarts<BoundariesSec(1)+60);%get the left licks that are within 1min of the trial start (that way you definitely get the licks after the stim has been delivered but you don't miss them if trial was terminated by manual reward delivery, for example)
    if isempty(LeftLick) %if no Licks, make it 0 to avoid problems of not being to index an empty array in if statements below
        LeftLick=0;
    else
        LeftLick=LeftLick-TrialStartTimes(trialNum); %convert to aligned to trial onset
        LeftLick=LeftLick'; %transpose (Eric's is a column)
    end
    Log{trialNum,12}={LeftLick};
    
    %Response
    if ismember('A', Stim) || (any(LeftLick==0) && any(RightLick==0)) %if there were premature licks and the stim wasn't delivered, the trial was aborted: 0
        Response=0;
    else
        FirstLick=min([RightLick(RightLick>(StimOnset) & RightLick<(StimOnset+1.5)); LeftLick(LeftLick>(StimOnset) & LeftLick<(StimOnset+1.5)) ]); % find the first lick within the response period
        if isempty(FirstLick)
            Response=0;
        elseif ismember(FirstLick, RightLick)
            Response=1;
        elseif ismember(FirstLick, LeftLick)
            Response=2;
        end
    end
    
    Log{trialNum,7}=Response;
    
    %Nex step is to incorporate the spikes/clusters
    
end

%% Trim and incorporate spiking data into Log

%Trim first and last 5 trials
Log=Log(6:end-5, :);
TrialStartTimes=TrialStartTimes(6:end-5);

%CRUCIAL: Remove the trials with no stims (premature licks for example)
temp=[]; %to build up the index of trials to keep
for i=1:length(Log)
    if ~isempty(Log{i,9}) %if there is NOT no stim_onset %%NOTE 20200210: this actually doesn't work because a cell with an empty array inside does not return 1 for isempty.. It doesn't hurt but is useless (the trials are actually dropped in the python script 'master_log'
        temp=[temp, i];
    end
end
Log=Log(temp, :); %keep only those trials
TrialStartTimes=TrialStartTimes(temp); %also update the start times array
%end of useless block

% get the spikes and incorporate into the Log
fnall = arrayfun(@(x) x.name(1:(end)), dir('*.t64'),'UniformOutput',false);

clusters = LoadSpikes(fnall);

mclustClusters=cell(1,length(clusters));
for i=1:length(clusters)
    mclustClusters{i}=data(clusters{i});
end

FinalLog=cell(length(TrialStartTimes)*length(mclustClusters), 14);

count=0;
for trial=1:size(Log,1)
    for unit=1:length(mclustClusters)
        count=count+1;
        
        %copy all info from Log
        FinalLog(count,:)=Log(trial,:);
        
        %fill column 13: spikes
        Aligned_spikes=mclustClusters{unit}-TrialStartTimes(trial); %align to trial start
        Aligned_trial_spikes=Aligned_spikes(Aligned_spikes>0 & Aligned_spikes<60);%to be very conservative, store 1min of spikes per trial
        FinalLog(count,13)={Aligned_trial_spikes};
        
        %fill column 14: unit name
        FinalLog(count,14)={['T', fnall{unit}(10:11), 'clst', fnall{unit}(14)]};
    end
end
log=FinalLog;%rename to be consistent with EF
save(['Log_', MouseName, '_', Date], 'log')

%Also add the structure a just to get the number of trials when
%calling PlotSingleUnitFromSingleLog_REVERSED.m
a=struct();
a.correctedIntanTrials=trial;
save('a', 'a')

%% Build Opto_log
optoTable=cell(length(clusters), 5);
optoTable(:,1)=log(1,1); %column1: MouseName
optoTable(:,2)=log(1,2); % column2: Date
for unit=1:length(clusters)
    optoTable(unit, 3)={['T', fnall{unit}(10:11), 'clst', fnall{unit}(14)]}; % column3: unit
end
optoTable(:, 4)={OptoStimStarts};
optoTable(:,5)={OptoStimEnds};

save(['Opto_log_', MouseName, '_', Date], 'optoTable')


%% Build waveform_log and optoSpikes_log

waveform_log = [];
optoSpikes_log = [];

for d = 1:size(optoTable,1)
    
    clustList = arrayfun(@(x) x.name(1:(end)), dir('*.t64'),'UniformOutput',false);
    StWv = Intan.clustered_ST_WV(clusters, clustList);
    
    endtime = TrialStartTimes(end);
    
    colNames = ['mouse_name', 'date', 'cluster_name', 'spikes', matlab.lang.makeUniqueStrings(repmat({'waveform'}, [1,128]))]; % (See note in RE64Create_Mouse_Waveforms_and_Timestamps.m)changed 128 to 256 because I now (starting ~7/2019) collect 64 timepoint per spike instead of 32 (so 32x4=128 replaced by 64x4=256
    
    cont_ind = cellfun(@(x) find(x<endtime, 150, 'last'), StWv(:,2), 'uni', 0);
    cont_clst = (cellfun(@(x, y,z) repmat({x}, [numel(z),1]), StWv(:,1), StWv(:,2), cont_ind, 'uni',0));
    cont_clst = vertcat(cont_clst{:});
    cont_st = cell2mat(cellfun(@(x, y) x(y), StWv(:,2), cont_ind, 'uni', 0));
    cont_wv = cell2mat(cellfun(@(x,y) x(y,:,:), StWv(:,3), cont_ind, 'uni', 0));
    
    if numel(cont_wv(1,1,:)) == 48
        cont_wv = [squeeze(cont_wv(:,1,8:39)), squeeze(cont_wv(:,2,8:39)),squeeze(cont_wv(:,3,8:39)),squeeze(cont_wv(:,4,8:39))];
    else
        cont_wv = [squeeze(cont_wv(:,1,:)), squeeze(cont_wv(:,2,:)),squeeze(cont_wv(:,3,:)),squeeze(cont_wv(:,4,:))];
    end
    
    cont_waves = [table(repmat(categorical({MouseName}),[numel(cont_clst),1]),...
        repmat(categorical({Date}), [numel(cont_clst),1]),...
        categorical(cont_clst), cont_st, 'VariableNames', colNames(1:4)),...
        array2table(single(cont_wv), 'VariableNames', colNames(5:end))];
    
    clst = cellfun(@(x, y, z) repmat({x}, [size(y(z(1):end),1),1]), StWv(:,1), StWv(:,2), cont_ind, 'uni',0);
    st = cellfun(@(x, y) x(y(1):end), StWv(:,2), cont_ind, 'uni', 0);
    wv = cell2mat(cellfun(@(x,y) x(y(1):end,:,:), StWv(:,3), cont_ind, 'uni', 0));
    if numel(wv(1,1,:)) == 48
        wv = [squeeze(wv(:,1,8:39)), squeeze(wv(:,2,8:39)),squeeze(wv(:,3,8:39)),squeeze(wv(:,4,8:39))];
    else
        wv = [squeeze(wv(:,1,:)), squeeze(wv(:,2,:)),squeeze(wv(:,3,:)),squeeze(wv(:,4,:))];
    end
    clst = vertcat(clst{:});
    st = vertcat(st{:});
    
    lightEvokedSp_ind = cell2mat(arrayfun(@(x) find(st>x & st<x+0.025), OptoStimStarts, 'uni',0)');
    %             lightEvokedSp_ind = lightEvokedSp_ind(unique(randi(numel(lightEvokedSp_ind), [150,1])));
    light_waves = [table(repmat(categorical({MouseName}),[numel(lightEvokedSp_ind),1]),...
        repmat(categorical({Date}), [numel(lightEvokedSp_ind),1]),...
        categorical(clst(lightEvokedSp_ind)), single(st(lightEvokedSp_ind)),...
        'VariableNames', colNames(1:4)),...
        array2table(single(wv(lightEvokedSp_ind,:)), 'VariableNames', colNames(5:end))];
    
    waves = [cont_waves; light_waves];
    optoSpikes = table(repmat(categorical({MouseName}),  [numel(st),1]), repmat(categorical({Date}),...
        [numel(st),1]), categorical(clst), st, 'VariableNames', colNames(1:4));
    
    clear StWv
    
    waveform_log = [waveform_log; waves];
    optoSpikes_log = [optoSpikes_log; optoSpikes];
end

    writetable(waveform_log, ['waveform_log_', MouseName, '_', Date, '.csv'])

    writetable(optoSpikes_log, ['optoSpikes_log_', MouseName, '_', Date, '.csv'])
%     waveform_log = table2cell(waveform_log);
     save(['waveform_log_', MouseName], 'waveform_log', '-v7')
%     save(['optoSpikes_log_', mouseName], 'optoSpikes_log')
%     clearvars -except workingDir mouseDir


%% ..........................................................................
%..........................................................................
% %%%%%%%%%% Below is only to plot behavior performance (partial copy )
% log=Log; %to plot behavior, use the Log (capital L) that has only 1 row per trial, no spike info
% TT=1;
% NumTT=1;
% %First, create intermediary indices that you can combine to pull out the
% %set of units or trials you want to plot
% idx_TT=[TT:NumTT:length(log(:,14))];
% %idx_Date=find(ismember(log(:,2),Date));
% idx_SomStim=find(ismember(log(:,4),'Stim_Som_NoCue'));
% idx_VisStim=find(ismember(log(:,4),'Stim_Vis_NoCue'));
% idx_NoStim=find(ismember(log(:,4), 'NAN')); %added 20190203 to discount the premature lick trials in calculation of performance
% idx_VisBlock=find(ismember(log(:,3),'Visual'));
% idx_VisBlock_Premature=intersect(idx_VisBlock, idx_NoStim);%added 20190203
% idx_SomBlock=find(ismember(log(:,3),'Whisker'));
% idx_SomBlock_Premature=intersect(idx_SomBlock, idx_NoStim);%added 20190203
% idx_Correct=find((ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==1) | ...
%     (ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==0) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==2) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==0));
% idx_False=find((ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==2) | ...
%     (ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==0) | ...
%     (ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==2) |...
%     (ismember(log(:,3),'Visual') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==1) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==1) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Som_NoCue') & cell2mat(log(:,7))==0) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==2) |...
%     (ismember(log(:,3),'Whisker') & ismember(log(:,4),'Stim_Vis_NoCue') & cell2mat(log(:,7))==1));
% 
% %idx_Unit=intersect(idx_TT, idx_Date); %This vector provides the index for all rows from the same unit
% idx_Unit=idx_TT; %In this SINGLE LOG version, no need to intersect with date
% idx_Unit_VisStim=intersect(idx_Unit, idx_VisStim); %This vector provides the index for all rows from the same unit dor Vis trials
% idx_Unit_SomStim=intersect(idx_Unit, idx_SomStim); %This vector provides the index for all rows from the same unit for SOM trials
% idx_Unit_VisBlock=intersect(idx_Unit, idx_VisBlock); %This vector provides the index for all rows from the same unit for Vis blocks
% idx_Unit_SomBlock=intersect(idx_Unit, idx_SomBlock); %This vector provides the index for all rows from the same unit for SOM blocks
% 
% idx_Unit_VisStim_VisBlock=intersect(idx_Unit_VisStim, idx_VisBlock); %This vector provides the index for all rows from the same unit for Vis trials on Vis blocks
% idx_Unit_VisStim_SomBlock=intersect(idx_Unit_VisStim, idx_SomBlock); %This vector provides the index for all rows from the same unit for Vis trials on Som blocks
% idx_Unit_SomStim_VisBlock=intersect(idx_Unit_SomStim, idx_VisBlock); %This vector provides the index for all rows from the same unit for SOM trials on Vis blocks
% idx_Unit_SomStim_SomBlock=intersect(idx_Unit_SomStim, idx_SomBlock); %This vector provides the index for all rows from the same unit for SOM trials on Som blocks
% 
% %% Now, arrange each by response type (by BlockType, StimType, Response)
% 
% %1st arrange them by Hit/Miss/FalseAlarm/CorrectRejection
% %Note: in the 'Response' column, 2 is lick left, 1 is lick right
% 
% % The purpose of the following block is to create two arrays. One with the
% % time of first lick after stim(regardless of side) along with the
% % trial/row number, and a second with the indices of NoResponse trials.
% 
% %For licks: the for loop below creates an array with 2 colums(firstlicktime,
% %rowinfullog) and as many rows as log
% LeftFirstLick=zeros(size(log, 1),2);
% RightFirstLick=zeros(size(log, 1),2);
% 
% 
% for k=1:size(log, 1)
%     if ~isempty(cell2mat([log{k,9}])) %This FOR/IF added on 20171217 because I ran into the problem of some trials having no StimOnsetTime (Column9)... not sure why these get incordporated from the Intan.sessionLog(a,b1) command but this IF statement should simply skip over them (needs to be tacked onto each loop that uses k (through the entire log)
%         
%         if isempty(cell2mat(log{k,12})) %if there was no lick, make it 0 such that it'll be sorted first
%             LeftFirstLick(k,1)=0;
%             LeftFirstLick(k,2)=k;
%         else
%             Licks=cell2mat(log{k,12})-cell2mat(log{k,9}); % if there are licks, substract stim time such that they are aligned with the spikes.
%             if isempty(Licks(Licks>0 & Licks<1.5)) %if all the licks are negtive, then there was no response -> 0
%                 LeftFirstLick(k,1)=0;
%                 LeftFirstLick(k,2)=k;
%             else
%                 LeftFirstLick(k,1)=min(Licks(Licks>0)); %Licks with negative values are before the stim, therefore don't count
%                 LeftFirstLick(k,2)=k;
%             end
%         end
%         
%         if isempty(cell2mat(log{k,11})) %if there was no lick, make it 0 such that it'll be sorted first
%             RightFirstLick(k,1)=0;
%             RightFirstLick(k,2)=k;
%             
%         else
%             Licks=cell2mat(log{k,11})-cell2mat(log{k,9});
%             if isempty(Licks(Licks>0 & Licks<1.5)) %if all the licks are negtive, then there was no response -> 0
%                 RightFirstLick(k,1)=0;
%                 RightFirstLick(k,2)=k;
%             else
%                 % if there are licks, substract stim time such that they are aligned with the spikes.
%                 RightFirstLick(k,1)=min(Licks(Licks>0)); %Licks with negative values are before the stim, therefore don't count
%                 RightFirstLick(k,2)=k;
%             end
%         end
%     end
% end
% 
% %The followinf loop goes through the possibilites and populates the
% %FirstLick ([time, idx]) and idx_NoResponse.
% FirstLick=zeros(size(log, 1),2);
% idx_NoResponse=[];
% for k=1:size(log, 1)
%     if ~isempty(cell2mat([log{k,9}])) %This FOR/IF added on 20171217 because I ran into the problem of some trials having no StimOnsetTime (Column9)... not sure why these get incordporated from the Intan.sessionLog(a,b1) command but this IF statement should simply skip over them
%         
%         if LeftFirstLick(k)~=0 || RightFirstLick(k)~=0
%             if LeftFirstLick(k)==0 && RightFirstLick(k)~=0
%                 Fir
stLick(k,1)=RightFirstLick(k);
%             elseif LeftFirstLick(k)~=0 && RightFirstLick(k)==0
%                 FirstLick(k,1)=LeftFirstLick(k);
%             elseif LeftFirstLick(k)> RightFirstLick(k)
%                 FirstLick(k,1)=RightFirstLick(k);
%             elseif LeftFirstLick(k)< RightFirstLick(k)
%                 FirstLick(k,1)=LeftFirstLick(k);
%             else
%                 FirstLick(k,1)=0;
%             end
%         else
%             idx_NoResponse=[idx_NoResponse;k];
%         end
%         FirstLick(k,2)=k;
%     end
% end
% TempLickIdx=FirstLick(:,1)~=0; %get a logical array identifying the FirstLick rows
% idx_FirstLick=FirstLick(TempLickIdx,2);% create the associated indices
% 
% 
% 
% % sum(FirstLick(:,1)~=0) %This was just a sanity check
% 
% %Need to create intermediate indices to categorize responses
% idx_Correct_Response=intersect(idx_Correct, idx_FirstLick);
% idx_Correct_NoResponse=intersect(idx_Correct, idx_NoResponse);
% idx_False_Response=intersect(idx_False, idx_FirstLick);
% idx_False_NoResponse=intersect(idx_False, idx_NoResponse);
% 
% %Now create the specific indices for each category
% idx_Unit_VisStim_VisBlock_Hit=intersect(idx_Unit_VisStim_VisBlock, idx_Correct_Response); %This vector provides the index for all rows from the same unit for Vis trials on Vis blocks for hits
% idx_Unit_VisStim_SomBlock_CorrectRejection=intersect(idx_Unit_VisStim_SomBlock, idx_Correct_NoResponse); %This vector provides the index for all rows from the same unit for Vis trials on Som blocks for correctrjections
% idx_Unit_VisStim_VisBlock_Miss=intersect(idx_Unit_VisStim_VisBlock, idx_False_NoResponse); %  no lick on VisStimVisBlock
% idx_Unit_SomBlock_FalseAlarm=intersect(idx_Unit_SomBlock, idx_False_Response); %This category lumps together 3 types of FalseAlamrs: rightORleft lick on VisStimSomBlock AND Right lick on VisStimVisBlock
% 
% idx_Unit_SomStim_SomBlock_Hit=intersect(idx_Unit_SomStim_SomBlock, idx_Correct_Response); %This vector provides the index for all rows from the same unit for Som trials on Som blocks for hits
% idx_Unit_SomStim_VisBlock_CorrectRejection=intersect(idx_Unit_SomStim_VisBlock, idx_Correct_NoResponse); %This vector provides the index for all rows from the same unit for Som trials on Vis blocks for correctrjections
% idx_Unit_SomStim_SomBlock_Miss=intersect(idx_Unit_SomStim_SomBlock, idx_False_NoResponse); %  no lick on SomStimSomBlock
% idx_Unit_VisBlock_FalseAlarm=intersect(idx_Unit_VisBlock, idx_False_Response); %This category lumps together 3 types of FalseAlamrs: rightORleft lick on SomStimVisBlock AND Left lick on SomStimSomBlock
% 
% 
% % PLOT
% hold on
% 
% scatter( idx_Unit_VisStim_VisBlock_Hit, ones(1,length( idx_Unit_VisStim_VisBlock_Hit)), '.g')
% scatter( idx_Unit_VisStim_SomBlock_CorrectRejection, ones(1,length( idx_Unit_VisStim_SomBlock_CorrectRejection)), '.k')
% scatter( idx_Unit_SomStim_SomBlock_Hit, ones(1,length( idx_Unit_SomStim_SomBlock_Hit))+1, '.g')
% scatter( idx_Unit_SomStim_SomBlock_Miss, ones(1,length( idx_Unit_SomStim_SomBlock_Miss))+1, '.k')
% scatter( idx_Unit_SomStim_VisBlock_CorrectRejection, ones(1,length( idx_Unit_SomStim_VisBlock_CorrectRejection))+1, '.k')
% scatter( idx_Unit_SomBlock_FalseAlarm, ones(1,length( idx_Unit_SomBlock_FalseAlarm)), '.r')
% scatter( idx_Unit_VisBlock_FalseAlarm, ones(1,length( idx_Unit_VisBlock_FalseAlarm))+1, '.r')
% scatter( idx_Unit_VisStim_VisBlock_Miss, ones(1,length( idx_Unit_VisStim_VisBlock_Miss)), '.k')
% 
% Touch=length(intersect(idx_Correct, idx_SomBlock))/(length(idx_SomBlock)-length(idx_SomBlock_Premature));
% Visual=length(intersect(idx_Correct, idx_VisBlock))/(length(idx_VisBlock)-length(idx_VisBlock_Premature));
% 
% title(['Correct Som: ', num2str(Touch), '    Correct Vis: ', num2str(Visual)])
% ylim([0, 4])
% 
% %save([MouseName, '_', Date, '_behavior.fig']);


