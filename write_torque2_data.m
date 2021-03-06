projpath = 'd:\Users\Gary\Google Drive\Muscle modeling\Min_jerk_files';
datafolder = [projpath filesep 'Data'];
cd(datafolder);
excel_file = 'sum_torque2_only.csv';
fileID=fopen(excel_file,'w');

delcount=0;

A={'c,'...
    'subj,'...
    's,'...
    't,'...
    'movedur,'...
    'sum_t2'...
    };

fprintf(fileID,'%s',A{:});
fprintf(fileID,'\n');

for c=1:4
    for subj = 1:8
        for s = 1:6
            for t = 1:8
                if ~isempty(shoulder{c,subj,s,t})
                    sum_t2 = (sum(shoulder{c,subj,s,t}.torque.^2)+sum(elbow{c,subj,s,t}.torque.^2))*...
                        vars{c,subj,s,t}.time_inc;
                    A={c,...
                        subj,...
                        s,...
                        t,...
                        vars{c,subj,s,t}.speeds(s),...
                        sum_t2,...
                        };
                        dlmwrite(excel_file,A,'-append')
                end
            end
        end
    end
end

fclose(fileID);