%script to run Responsiveness.m in multiple folders in a row
fMouse={'Claustrum18','Claustrum23','Claustrum25'};
fMouse={'Claustrum4','Claustrum5', 'Claustrum6'};
fMouse={'Claustrum31','Claustrum32', 'Claustrum37'};

for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    cd(MouseName)
    fDate = arrayfun(@(x) x.name(1:(end)), dir('*-17'),'UniformOutput',false);
    for i=1:length(fDate)
        cd(fDate{i})
        LogFile=arrayfun(@(x) x.name(1:(end)), dir('Log_*'), 'UniformOutput', false);
        OptoLog=arrayfun(@(x) x.name(1:(end)), dir('Opto_log_*'), 'UniformOutput', false);
        Range_B=0.4;
        Salt_window=0.01;
        Mouse=MouseName;
        Freq_Block='2Hz';
        Output_location='C:\Users\Brown Lab\Downloads\';
        OptoNetworkTagging2_SALT(Mouse,LogFile{1},OptoLog{1}, Range_B, Salt_window,Freq_Block, Output_location )
        disp(fDate{i})
        cd ..
    end
    cd ..
end
