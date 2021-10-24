
p_values=[];
probas=[];
%script to run Responsiveness.m in multiple folders in a row
fMouse={'Claustrum4','Claustrum5','Claustrum5'};
for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    cd(MouseName)
    fDate = arrayfun(@(x) x.name(1:(end)), dir('*-17'),'UniformOutput',false);
    for i=1:length(fDate)
        cd(fDate{i})
        LogFile=arrayfun(@(x) x.name(1:(end)), dir('Log_*'), 'UniformOutput', false);
        a='a.mat';%same for everyone
        OptoLog=arrayfun(@(x) x.name(1:(end)), dir('Opto_log_*'), 'UniformOutput', false);
        Range_B=0.4;
        Salt_window=0.02
        %Note on SALT: the Range influences the resolution of you p-value.
        %If you pick a range of 50msec and a test window of 10msec
        %(default), the resolution of your p-value will be 0.1 (5 null
        %windows for 1 test window, 4+3+2+1=10 pairwise distances can be
        %calculated so the probability of your test compared to this
        %distribution is 0, 0.1, 0.2, .... So to increase the resolution,
        %increase the length of the baseline window.
        Mouse=MouseName;
        [p,prob]=OptoNetworkTagging2_SALT_savepvalues(Mouse,LogFile{1}, a,OptoLog{1}, Range_B, Salt_window);
        p_values=[p_values,p];
        probas=[probas,prob];
        disp(fDate{i})
        cd ..
    end
    cd ..
end

%script to run Responsiveness.m in multiple folders in a row
fMouse={'Claustrum31','Claustrum32','Claustrum37'};
for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    cd(MouseName)
    fDate = arrayfun(@(x) x.name(1:(end)), dir('2019*'),'UniformOutput',false);
    for i=1:length(fDate)
        cd(fDate{i})
        LogFile=arrayfun(@(x) x.name(1:(end)), dir('Log_*'), 'UniformOutput', false);
        a='a.mat';%same for everyone
        OptoLog=arrayfun(@(x) x.name(1:(end)), dir('Opto_log_*'), 'UniformOutput', false);
        Range_B=0.4;
        Salt_window=0.02;
        Mouse=MouseName;
        [p,prob]=OptoNetworkTagging2_SALT_savepvalues(Mouse,LogFile{1}, a,OptoLog{1}, Range_B, Salt_window);
        p_values=[p_values,p];
        probas=[probas,prob];
        disp(fDate{i})
        cd ..
    end
    cd ..
end


%script to run Responsiveness.m in multiple folders in a row
fMouse={'Claustrum18','Claustrum23','Claustrum25'};
for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    cd(MouseName)
    fDate = arrayfun(@(x) x.name(1:(end)), dir('2019*'),'UniformOutput',false);
    for i=1:length(fDate)
        cd(fDate{i})
        LogFile=arrayfun(@(x) x.name(1:(end)), dir('Log_*'), 'UniformOutput', false);
        a='a.mat';%same for everyone
        OptoLog=arrayfun(@(x) x.name(1:(end)), dir('Opto_log_*'), 'UniformOutput', false);
        Range_B=0.4;
        Salt_window=0.02;
        Mouse=MouseName;
        [p,prob]=OptoNetworkTagging2_SALT_savepvalues(Mouse,LogFile{1}, a,OptoLog{1}, Range_B, Salt_window);
        p_values=[p_values,p];
        probas=[probas,prob];
        disp(fDate{i})
        cd ..
    end
    cd ..
end

%Figure
fill(log([0.0019,0.05, 0.05, 0.0019]),log([0.1,0.1,2,2]),'r', 'LineStyle','none')
alpha(0.2)
hold on
scatter(log(p_values), log(probas),'b.')
for i=1:length(p_values)
    if p_values(i)==0
        scatter(log(p_values(i)+0.002), log(probas(i)),'b.')
    end
end
xlabel('SALT p-values')
xticks(log([0.05, 0.1,0.2,0.4]))
xticklabels({'0.05','0.1','0.2','0.4'})
ylabel('Mean number of spikes within 20msec of laser pulse')
yticks(log([0.1, 0.2,0.4,0.8, 1.6]))
yticklabels({'0.1','0.2','0.4', '0.8','1.6'})
vline(log(0.05))
hline(log(0.1))
xlim([log(0.0019),log(1)])

%counting..
count=0;
for i=1:length(p_values)
    if p_values(i)<0.05 && probas(i)>0.1
        count=count+1;
    end
end

