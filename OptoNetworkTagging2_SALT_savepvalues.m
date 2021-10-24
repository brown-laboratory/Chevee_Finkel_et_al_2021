function  [p_values, probas] = OptoNetworkTagging2_SALT_savepvalues(Mouse,LogFile, a,OptoFile, Range_B, Salt_window)
% This function performs SALT test to determine whether the distribution of
% first spike latencies is different pre/post laser stim.
%
% Based on Opto_SALT_test.m (calls salt.m, from Adam Kepecs' lab in CSH)
%
%Range (seconds): the window pre and post to test for a difference in spike counts
%(eg. 0.01)
p_values=[];
probas=[];
load(LogFile) %loading the file should create the variable log, which is a cell array with everthing...
%load(a);
load(OptoFile);
optoTable{1,6}=[]; %create the cells to fill in with the spike times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get spikes from time range during OptoStim protocol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnall = arrayfun(@(x) x.name(1:(end)), dir('*.t64'),'UniformOutput',false);
clusters = LoadSpikes(fnall);

mclustClusters=cell(1,length(clusters));
for i=1:length(clusters)
    mclustClusters{i}=data(clusters{i});
end

for i=1:length(mclustClusters)
    optoTable{i,6}=[mclustClusters{i}(mclustClusters{i}>=optoTable{i,4}(1))]'; %get spikes after the optostim protocol started
end


%Define Blocks
OptoStarts=optoTable(1,4); % create cell array with onset times
OptoEnds=optoTable(1,5);
% next identify the blocks
i=1; %initialize counter for identofying blocks
BlockStarts=[];
for i=3: length(OptoStarts{1})-1
    if (round([OptoStarts{1}(i)-OptoStarts{1}(i-1)],3)~=round([OptoStarts{1}(i+1)-OptoStarts{1}(i)], 3)) && (round([OptoStarts{1}(i+1)-OptoStarts{1}(i)],3)~=round([OptoStarts{1}(i-1)-OptoStarts{1}(i-2)],3)) && (round([OptoStarts{1}(i)-OptoStarts{1}(i-1)],3)~=round([OptoStarts{1}(i-1)-OptoStarts{1}(i-2)],3)) %This strategy finds the block starts that if the delay is different from both before and after it mu
        BlockStarts=[BlockStarts, i]; % this array now contains the indices to the optostim onsets
    end
end


for TT=1:length(mclustClusters)
    
    SPT_TEST=zeros(length(BlockStarts(1):BlockStarts(2)-1),Range_B*30000); %this creates the backbone for the intput test , we are using a test window of 50msec
    %matrix to SALT. Number of rows is equal to number f=of pulses, number of
    %columsn is equal to the time wuindown in binary (time in secondxSampling raate)
    SPT_BASELINE=SPT_TEST; % ame backbone, but represents the window before the pulses
    
    %Because we are using a baseline of 200msec, we cannot use the 10Hz
    %pules, so we use the 2Hz
    j=0; %make the test matrix
    for i= optoTable{TT,4}(BlockStarts(1):BlockStarts(2)-1)
        j=j+1;
        temp=optoTable{TT,6}(optoTable{TT,6}>i & optoTable{TT,6}<i+Range_B)-i; %takes all the spikes in the 50ms following each pulse
        for m=ceil(temp*30000) %take each spike, convert to its binary position based on sampling rate (30000) and make a "1" in SPT matrix
            SPT_TEST(j,m)=1;
        end
    end
    
    
    j=0; %make the baseline matrix
    for i= optoTable{TT,4}(BlockStarts(1):BlockStarts(2)-1)
        j=j+1;
        temp=optoTable{TT,6}(optoTable{TT,6}<i & optoTable{TT,6}>i-Range_B)-(i-Range_B); %takes all the spikes in the 50ms following each pulse
        for m=ceil(temp*30000) %take each spike, convert to its binary position based on sampling rate (30000) and make a "1" in SPT matrix
            SPT_BASELINE(j,m)=1;
        end
    end
    
    
    [p,~]=salt(SPT_BASELINE, SPT_TEST, 1/30000, Salt_window);
    
    % In addition to SALT, we need to calculate the probablity and set a
    % threshold (as with for hard core optotag) because SALT is very sensitive
    % and will be triggered by a single spike.
    Proba=sum(sum(SPT_TEST(:,1:Salt_window*30000)))/length((BlockStarts(1):BlockStarts(2)-1)); %represents the number of spikes per 50msec post laser
    
    %Mouse=optoTable{TT,1};
    Date=optoTable{TT,2};
    Tetrodeclst=optoTable{TT,3};
    
    disp( Tetrodeclst)
    p_values=[p_values,p];
    probas=[probas,Proba];
%      if p<0.05 && Proba>0.1%ie. if the distribution of spikes is different and there are enough spikes
%         disp( Tetrodeclst)
%         fid = fopen(['E:\work form home 20200313\', 'test.txt'],'a');
%         %fid = fopen(['Z://Maxime Chevee/Recording DATA/',Mouse(1:9), ' ',Mouse(10:end), '/', Mouse,'_', 'OptoNetworkList_SALT.txt'],'a');
%         fprintf(fid, [Mouse '_' Date '_' Tetrodeclst '\n']);
%         fclose(fid);
%     end
end
