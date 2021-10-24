

fMouse={'Claustrum18','Claustrum23','Claustrum25','Claustrum31','Claustrum32','Claustrum37'};
fMouse={'Claustrum4','Claustrum5','Claustrum6'};

for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    cd(MouseName)
    f_target = fopen([MouseName,'_OptoNetworkList_SALT.txt'],'rt');
    t_line = fgetl(f_target);
    while ischar(t_line)
        f_rest= fopen([MouseName,'_RestList.txt'],'rt');
        r_line=fgetl(f_rest);
        while ischar(r_line)
            if all(t_line(end-17:end-5)==r_line(end-17:end-5))  %all(t_line(1:end-5)==r_line(1:end-5))
                if ~all(t_line(end-17:end)==r_line(end-17:end))  %~all(t_line==r_line)
                    disp(t_line)
                    fid = fopen(['Z://Maxime Chevee/Recording DATA/',MouseName, '/', MouseName,'_', 'OptoNetworkList_SALT_SameTT.txt'],'a');
                    fprintf(fid, [r_line, '\n']);
                    fclose(fid);
                end
            end
            r_line=fgetl(f_rest);
        end
        fclose(f_rest);
        t_line = fgetl(f_target);
    end
    fclose(f_target);  
    cd ..
end




% fixing the lists based on new SALT p-values (2020-04-01)


fMouse={'Claustrum18','Claustrum23','Claustrum25','Claustrum31','Claustrum32','Claustrum37'};
fMouse={'Claustrum4','Claustrum5','Claustrum6'};

for mouse=1:length(fMouse)
    MouseName=fMouse{mouse};
    f_target = fopen([MouseName,'_OptoNetworkList_SALT.txt'],'rt');
    t_line = fgetl(f_target);
    while ischar(t_line)
        f_rest= fopen([MouseName,'_RestList.txt'],'rt');
        r_line=fgetl(f_rest);
        while ischar(r_line)
            if  all(t_line(1:end-5)==r_line(1:end-5))%all(t_line(end-17:end-5)==r_line(end-17:end-5))%   %
                if ~all(t_line==r_line)%~all(t_line(end-17:end)==r_line(end-17:end))%   % 
                    disp(t_line)
                    fid = fopen(['E:\work form home 20200313\', MouseName,'_', 'OptoNetworkList_SALT_SameTT.txt'],'a');
                    fprintf(fid, [r_line, '\n']);
                    fclose(fid);
                end
            end
            r_line=fgetl(f_rest);
        end
        fclose(f_rest);
        t_line = fgetl(f_target);
    end
    fclose(f_target);  
end



