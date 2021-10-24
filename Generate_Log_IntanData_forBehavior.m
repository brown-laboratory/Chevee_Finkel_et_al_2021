function Generate_Log_IntanData_forBehavior(MouseName, Date)

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
save(['Log_',MouseName, '_', Date, '.mat'], 'Log');
save(['IntanData_',MouseName, '_', Date, '.mat'], 'RightRewardStarts', 'LeftRewardStarts', 'RightLickStarts', 'LeftLickStarts', 'WhiskerStimStarts', 'VisualStimStarts', 'VisualBlockStarts','WhiskerBlockStarts');
