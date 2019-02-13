%Load Undorded without FML

projpath = 'D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\'
datafolder = [projpath filesep 'Data'];
cd(datafolder);
excel_file = 'sum_torque2_only.csv';
fileID=fopen(excel_file,'w');

delcount=0;

A={'c,'...
    's,'...
    't,'...
    'movedur,'...
    'sum_t2'...
    };

fprintf(fileID,'%s',A{:});
fprintf(fileID,'\n');

for t = 1:8
    for c=1:4
        for s = 1:6
            sum_t2 = sum(shoulder{c,s,t}.torque.^2)+sum(elbow{c,s,t}.torque.^2);
            [c,s,sum_t2]
            A={c,...
                s,...
                t,...
                vars{c,s,t}.speeds(s),...
                sum_t2,...
                };
                dlmwrite(excel_file,A,'-append')             
        end
    end
end

fclose(fileID);

close all
close all
close all
close all
close all
close all
close all