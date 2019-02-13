minparams = {'umberger','stress','force','act','drive','stress2','force2','act2','drive2'};
% minparams = {'stress','force','act','drive','stress2','force2','act2','drive2'};

main_folder = 'D:\Users\Gary\Desktop\Model Testing';
count = 1;
for k1 = 1:length(minparams)
    cd([main_folder filesep 'Data']);
    load(sprintf('aa_%scost_05-07-2018.mat',minparams{k1}));
    for i =1:4*8*7*8

        [c,subj,s,t] = ind2sub([4,8,7,8],i);

        shoulder1{c,subj,s,t} = shoulder{i};
        elbow1{c,subj,s,t} = elbow{i};
        theta1{c,subj,s,t} = theta{i};
        eff_mass1{c,subj,s,t} = eff_mass{i};
        muscles1{c,subj,s,t} = muscles{i};
        act1{c,subj,s,t} = act{i};
        u1{c,subj,s,t} = u{i};
        est1{c,subj,s,t} = est{i};
        tnew1{c,subj,s,t} = tnew{i};
        energy1{c,subj,s,t} = energy{i};

        shoulder{i} = [];
        elbow{i} = [];
        theta{i} = [];
        eff_mass{i} = []; 
        muscles{i} = [];
        act{i} = [];
        u{i} = [];
        est{i} = [];
        tnew{i} = [];
        energy{i} = [];

    end
    clear shoulder elbow theta eff_mass muscles act u est tnew energy
    for i =1:4*8*7*8

        [c,subj,s,t] = ind2sub([4,8,7,8],i);

        shoulder{c,subj,s,t} = shoulder1{c,subj,s,t};
        elbow{c,subj,s,t} = elbow1{c,subj,s,t};
        theta{c,subj,s,t} = theta1{c,subj,s,t};
        eff_mass{c,subj,s,t} = eff_mass1{c,subj,s,t};
        muscles{c,subj,s,t} = muscles1{c,subj,s,t};
        act{c,subj,s,t} = act1{c,subj,s,t};
        u{c,subj,s,t} = u1{c,subj,s,t};
        est{c,subj,s,t} = est1{c,subj,s,t};
        tnew{c,subj,s,t} = tnew1{c,subj,s,t};
        energy{c,subj,s,t} = energy1{c,subj,s,t};

        shoulder1{c,subj,s,t} = [];
        elbow1{c,subj,s,t} = [];
        theta1{c,subj,s,t} = [];
        eff_mass1{c,subj,s,t} = []; 
        muscles1{c,subj,s,t} = [];
        act1{c,subj,s,t} = [];
        u1{c,subj,s,t} = [];
        est1{c,subj,s,t} = [];
        tnew1{c,subj,s,t} = [];
        energy1{c,subj,s,t} = [];

    end
    filename = sprintf('aa_%scost_05-07-2018',  minparams{L});
    % if exist('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Data')==7
    %     cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Data');
    % else
    %     cd('C:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Data');    
    % end
    cd([fold '\Data']);
    save(filename);
    if exist(strcat(filename,'.mat'),'file')==2
        fprintf('File saved as %s \n',filename');
    end 
end