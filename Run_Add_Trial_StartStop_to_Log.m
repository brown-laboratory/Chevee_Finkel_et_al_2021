 %script to runAdd_Trial_StartStop_to_Log(MouseName, Date) in multiple folders in a row
 
 fMouse={'Claustrum37'};
 for mouse=1:length(fMouse)
     MouseName=fMouse{mouse};
     cd(MouseName)
     MouseName=[MouseName(1:9), MouseName(end-1:end)]; % this ensures tat if the folder has a space between 'Claustrum' and '37' it will be dropped when naming 'Log'
     fDate = arrayfun(@(x) x.name(1:(end)), dir('2019*'),'UniformOutput',false);
     for i=1:length(fDate)
         Date=fDate{i};
         cd(Date)
         
         Add_Trial_StartStop_to_Log(MouseName, Date)
         disp(fDate{i})
         cd ..
     end
     cd ..
 end
 
 

