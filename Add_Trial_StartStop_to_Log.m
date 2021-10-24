function Add_Trial_StartStop_to_Log(MouseName, Date)
%This function serves to add a 15th column to the Log.
% It is only useful for mice processed on brown lab rig, because in EF
% version of log the Licks are limited to within the trial, in mine they
% are fom start to 60sec, which means I don't have the info of when a trial
% stops. I will add this 15th column to the script Generate_log_EFTask_BrownLabRig_MClust(MouseName, Date)
% and EFtask_analysis_BrownLabRig_MClust.m so that this script shouldn't be
% needed again.

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
            aux_input_channels,aux_input_data,t_aux_input, t_dig, board_dig_in_data, board_adc_data]= Intan.read_Intan_RHD2000_file_noUI(file(2:end-2));%NOTE: depending on MATLAB version, the indexing is slightly off, try file(3:end-3) 
        numChan=numel(amplifier_channels);
        fs=frequency_parameters.amplifier_sample_rate;
    else
        [notes,frequency_parameters,amplifier_channels,amplifier_data,t_amplifier,...
            aux_input_channels,aux_input_data,t_aux_input, t_dig, board_dig_in_data, board_adc_data]= Intan.read_Intan_RHD2000_file_noUI(file(2:end-2));%NOTE: depending on MATLAB version, the indexing is slightly off, try file(3:end-3) 
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
Log=cell(length(TrialStartTimes), 1);

%these are the temporary attributes of the current block, we will acquire
%them and build the row at the end

trialNum=0;
%for i=1:length(board_dig_in-data); % we will loop through each 1/30000th of a second and update the current conditions as we go
for i=1:length(TrialStartTimes)-1
    trialNum=i;
    %get the boudaries of the trial
    BoundariesSec=[TrialStartTimes(trialNum), TrialStartTimes(trialNum+1)]; %start and stop of trial in seconds
Log{trialNum}=BoundariesSec;
end
Log=Log(6:end-5, :);
TTs=arrayfun(@(x) x.name(1:(end)), dir('*.t64'), 'UniformOutput', false);
LogFile=arrayfun(@(x) x.name(1:(end)), dir('Log_Claustrum*'), 'UniformOutput', false);
x=load(LogFile{1}); %load the original 'log'
log=x.log;
count=0;
for row=1:length(TTs):length(log)
    count=count+1;
        log(row:row+length(TTs)-1,15)=deal(Log(count));
end
save(['Log_', MouseName, '_', Date], 'log')


