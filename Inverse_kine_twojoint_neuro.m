         
 %% Init vars
% 
% cd('D:\Users\Gary\OneDrive\Muscle modeling\Metabolics Files');
% if exist('Data')~=1
%     load('XY_outstop_data.mat');
%     fprintf('Loaded data \n');
% end
clc
%% Init Vars
% clear
% global time_step
time_step = .0025;
% global minfail

target_str = {'one','two','thr','fou','fiv','six','sev','eig'};
input_normforce = 200E4;%31.84E4;%200E4
% masses = [0,5,10,20]/2.2;
minparams = {'umberger','stress','force','act','drive','stress2','force2','act2','drive2'};
% minparams = {'umberger','act','drive','stress2','force2','act2','drive2'};

% global t_act ;
t_act= 0.05;%.0050;
% global t_deact;
t_deact = 0.066;%.066;

fold = pwd;
% cd([fold '\Data']);
cd('D:\Users\Gary\Desktop\Model Testing\Data');
load('Resamp_data_gb.mat');
cd(fold);

for L = 1:length(minparams)
tic
minfail = [];
filename = sprintf('aa_%scost_05-07-2018',  minparams{L});

fprintf('Processing %s \n',filename);
clearvars -except time_step target_str input_normforce masses minparams L Resamp tic fold minfail
subj_mass = [125,156.199,194.7,129.8,173.8,159.7,136.4,127.6];
subj_heights = [162.5,167.60,195.6,167.6,180,172.7,175,171.5]/100;
thetaE_start = [1.7049,1.7588,1.9786,1.7588,1.8696,1.8075,1.8280,1.7965];
thetaS_start = [0.7008,0.6644,0.5089,0.6644,0.5876,0.6309,0.6167,0.6386];

upperarm_mass = 0.028*subj_mass/2.2;
lowerarm_mass = 0.022*subj_mass/2.2;
upperarm_length = (.825-.632)*subj_heights;
lowerarm_length = (.632-0.370)*subj_heights;
upperarm_lcom = 0.436*upperarm_length;
lowerarm_lcom = 0.682*lowerarm_length;

for c=1:4
    for subj=1:8
    for s=1:7
        for t = 1:8
            for d = 2:2
    if sum(isnan(Resamp.P(c,subj,s,:)))>10
        continue
    end
    
    vars{c,subj,s,t}.time_inc = time_step;
    vars{c,subj,s,t}.speeds = squeeze(Resamp.Endt(c,subj,:,t));
    switch c
        case 1
            vars{c,subj,s,t}.masses = 0/2.2;
        case 2
            vars{c,subj,s,t}.masses = 5/2.2;
        case 3
            vars{c,subj,s,t}.masses = 10/2.2+1;
        case 4
            vars{c,subj,s,t}.masses = 20/2.2;
    end
    
    vars{c,subj,s,t}.distances = [.05,.1,.15,.2];
    vars{c,subj,s,t}.norm_force = input_normforce;%200E4;%31.8E4;
    vars{c,subj,s,t}.minparam = minparams{L};
    
    forearm{c,subj,s,t}.mass = lowerarm_mass(subj);%;+2*vars{c,subj,s,t}.masses(c);
    forearm{c,subj,s,t}.length = lowerarm_length(subj); % meters, taken from An iterative optimal control and estimation design for nonlinear stochastic system
    forearm{c,subj,s,t}.l_com = 0.682*lowerarm_length(subj);
    [forearm{c,subj,s,t}] = calc_forearmI(forearm{c,subj,s,t},vars{c,subj,s,t}.masses);
    
    upperarm{c,subj,s,t}.length = upperarm_length(subj); % meters
    upperarm{c,subj,s,t}.l_com = 0.436*upperarm_length(subj);
    upperarm{c,subj,s,t}.centl = 0.436*upperarm_length(subj);
    upperarm{c,subj,s,t}.mass = upperarm_mass(subj); %kg
    upperarm{c,subj,s,t}.Ic = .0141;

    shoulder{c,subj,s,t} = [];
    elbow{c,subj,s,t} = [];
    theta{c,subj,s,t} = [];
    
%     ro = [0,.4];
    center = [-.0758,0.4878];
    switch t
        case 1
            vars{c,subj,s,t}.targets.one = [.0707 .0707];
            ro = [-.0758,0.4878];
        case 2
            vars{c,subj,s,t}.targets.two = [-.0707 .0707];
            ro = [-.0758,0.4878];
        case 3
            vars{c,subj,s,t}.targets.thr = [-.0707 -.0707];
            ro = [-.0758,0.4878];
        case 4
            vars{c,subj,s,t}.targets.fou = [.0707 -.0707];
            ro = [-.0758,0.4878];
        case 5
            ro = [.0707 .0707]+center;
            vars{c,subj,s,t}.targets.fiv = [-.0707 -.0707];
        case 6
            ro = [-.0707 .0707]+center;
            vars{c,subj,s,t}.targets.six = [.0707 -.0707];
        case 7
            ro = [-.0707 -.0707]+center;
            vars{c,subj,s,t}.targets.sev = [.0707 .0707];
        case 8
            ro = [.0707 -.0707]+center;
            vars{c,subj,s,t}.targets.eig = [-.0707 .0707];
    end
            
    
    
    rf = [vars{c,subj,s,t}.targets.(target_str{t})(1)*vars{c,subj,s,t}.distances(d)/.1+ro(1),...
        vars{c,subj,s,t}.targets.(target_str{t})(2)*vars{c,subj,s,t}.distances(d)/.1+ro(2)];
    
%     [Data{c,subj,s,t}.time,Data{c,subj,s,t}.x,Data{c,subj,s,t}.y,~,...
%         ~,~,~] =...
%         Gen_mvt_gb(Resamp,c,s,t,ro,rf,vars{c,subj,s,t}.time_inc);
    
    [Data{c,subj,s,t}] = Gen_mvt_gb(Resamp,c,subj,s,t,ro,rf,vars{c,subj,s,t}.time_inc);
    
%     [Data{c,subj,s,t}.time,Data{c,subj,s,t}.x,Data{c,subj,s,t}.y,Data{c,subj,s,t}.vx,...
%         Data{c,subj,s,t}.vy,Data{c,subj,s,t}.ax,Data{c,subj,s,t}.ay] =...
%         Gen_mvt_gb(Resamp,c,subj,s,t,ro,rf,vars{c,subj,s,t}.time_inc);
    
    
%     [Data{c,subj,s,t}.time,Data{c,subj,s,t}.x,Data{c,subj,s,t}.y,Data{c,subj,s,t}.vx,...
%         Data{c,subj,s,t}.vy,Data{c,subj,s,t}.ax,Data{c,subj,s,t}.ay] ...
%         =minjerk(ro,rf,vars{c,subj,s,t}.speeds(s),vars{c,subj,s,t}.time_inc);
    
    Data{c,subj,s,t}.targetposition(1)=vars{c,subj,s,t}.targets.(target_str{t})(1)*vars{c,subj,s,t}.distances(d)/.1+ro(1);
    Data{c,subj,s,t}.targetposition(2)=vars{c,subj,s,t}.targets.(target_str{t})(2)*vars{c,subj,s,t}.distances(d)/.1+ro(2);
    
    Data{c,subj,s,t}.startposition = ro;
    
    vars{c,subj,s,t}.masses;
            end
        end
    end
    end
end
clear c subj s t
parfor_progress(4*8*7*8);
% tic
for i=1:4*8*7*8
    tic
    [c,subj,s,t] = ind2sub([4,8,7,8],i);
%     i
    if isempty(Data{c,subj,s,t})
%         parfor_progress;
        continue
    end
    [ shoulder1{i},elbow1{i},theta1{i},muscles1{i},act1{i},u1{i},est1{i},tnew1{i}, energy1{i}, eff_mass1{i} ] =...
        looper2(Data{c,subj,s,t}, forearm{c,subj,s,t}, upperarm{c,subj,s,t},vars{c,subj,s,t},time_step,c,subj,s,t);
%     [ shoulder1{i},elbow1{i},theta1{i}, eff_mass1{i}] =...
%         looper3(Data{c,subj,s,t}, forearm{c,subj,s,t}, upperarm{c,subj,s,t},vars{c,subj,s,t},time_step,c,s,t);
    parfor_progress;
%     disp(i)
toc;
1;
end
parfor_progress(0);

for i =1:4*8*7*8
   
    [c,subj,s,t] = ind2sub([4,8,7,8],i);
    
    shoulder{c,subj,s,t} = shoulder1{i};
    elbow{c,subj,s,t} = elbow1{i};
    theta{c,subj,s,t} = theta1{i};
    eff_mass{c,subj,s,t} = eff_mass1{i};
    muscles{c,subj,s,t} = muscles1{i};
    act{c,subj,s,t} = act1{i};
    u{c,subj,s,t} = u1{i};
    est{c,subj,s,t} = est1{i};
    tnew{c,subj,s,t} = tnew1{i};
    energy{c,subj,s,t} = energy1{i};
    
    shoulder1{i} = [];
    elbow1{i} = [];
    theta1{i} = [];
    eff_mass1{i} = []; 
    muscles1{i} = [];
    act1{i} = [];
    u1{i} = [];
    est1{i} = [];
    tnew1{i} = [];
    energy1{i} = [];
    
end
clear shoulder1 elbow1 theta1 muscles1 act1 est1 tnew1

% if ~isempty(minfail)
%     fprintf('The sim failed at: c,subj,s,t\n');
%     minfail
% end

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
cd(fold);
toc 
end
