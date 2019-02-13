         
 %% Init vars
% 
% cd('D:\Users\Gary\OneDrive\Muscle modeling\Metabolics Files');
% if exist('Data')~=1
%     load('XY_outstop_data.mat');
%     fprintf('Loaded data \n');
% end

%% Init Vars
% clear
% global time_step
time_step = .0025;

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
cd([fold '\Data']);
% cd('E:\Documents\Google Drive\Muscle modeling\Min_jerk files\Data');
load('Resamp_data_aa.mat');
cd(fold);

for L = 7:9%:length(minparams)
tic
filename = sprintf('aa_%scost_03-22-2018',  minparams{L});

fprintf('Processing %s \n',filename);
clearvars -except time_step target_str input_normforce masses minparams L Resamp tic fold
for c=1:4
    for s=1:6
        for t = 1:8
            for d = 2:2
                
    vars{c,s,t}.time_inc = time_step;
    switch c
        case 1
%         vars{c,s,t}.speeds = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
%         vars{c,s,t}.speeds = [0.4660, 0.5364, 0.6401, 0.8162, 1.0024, 1.1934];
        vars{c,s,t}.speeds = [0.3285, 0.4022, 0.5363, 0.7509, 0.9853, 1.1957];
        vars{c,s,t}.masses = 0/2.2;
        case 2
%         vars{c,s,t}.speeds = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
%         vars{c,s,t}.speeds = [0.4875, 0.5475, 0.6379, 0.8156, 0.9975, 1.1913];
        vars{c,s,t}.speeds = [0.3567, 0.4252, 0.5369, 0.7539, 0.9693, 1.1938];
        vars{c,s,t}.masses = 5/2.2;
        case 3
%         vars{c,s,t}.speeds = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
%         vars{c,s,t}.speeds = [.5473    0.6418    0.8162    0.9981     1.1823    1.3705];
        vars{c,s,t}.speeds = [0.4244, 0.5449, 0.7561, 0.9750, 1.1892, 1.3996];
        vars{c,s,t}.masses = 10/2.2+1;
        case 4
%         vars{c,s,t}.speeds = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
%         vars{c,s,t}.speeds = [0.5626    0.6440    0.8178    0.9960     1.1832    1.3584];
        vars{c,s,t}.speeds = [0.4372, 0.5481, 0.7647, 0.9723, 1.1915, 1.3896];
        vars{c,s,t}.masses = 20/2.2;
    end
    vars{c,s,t}.distances = [.05,.1,.15,.2];
    vars{c,s,t}.norm_force = input_normforce;%200E4;%31.8E4;
    vars{c,s,t}.minparam = minparams{L};
    
    forearm{c,s,t}.mass = 1.5186;%;+2*vars{c,s,t}.masses(c);
    forearm{c,s,t}.length = .33*1.3; % meters, taken from An iterative optimal control and estimation design for nonlinear stochastic system
    forearm{c,s,t}.l_com = (.33*1.3)*2/3;
    [forearm{c,s,t}] = calc_forearmI(forearm{c,s,t},vars{c,s,t}.masses);
    
    upperarm{c,s,t}.length = .33; % meters
    upperarm{c,s,t}.l_com = .33/2;
    upperarm{c,s,t}.centl = .33/2;
    upperarm{c,s,t}.mass = 1.93284; %kg
    upperarm{c,s,t}.Ic = .0141;

    shoulder{c,s,t} = [];
    elbow{c,s,t} = [];
    theta{c,s,t} = [];
    
%     ro = [0,.4];
    center = [-.0758,0.4878];
    switch t
        case 1
            vars{c,s,t}.targets.one = [.0707 .0707];
            ro = [-.0758,0.4878];
        case 2
            vars{c,s,t}.targets.two = [-.0707 .0707];
            ro = [-.0758,0.4878];
        case 3
            vars{c,s,t}.targets.thr = [-.0707 -.0707];
            ro = [-.0758,0.4878];
        case 4
            vars{c,s,t}.targets.fou = [.0707 -.0707];
            ro = [-.0758,0.4878];
        case 5
            ro = [.0707 .0707]+center;
            vars{c,s,t}.targets.fiv = [-.0707 -.0707];
        case 6
            ro = [-.0707 .0707]+center;
            vars{c,s,t}.targets.six = [.0707 -.0707];
        case 7
            ro = [-.0707 -.0707]+center;
            vars{c,s,t}.targets.sev = [.0707 .0707];
        case 8
            ro = [.0707 -.0707]+center;
            vars{c,s,t}.targets.eig = [-.0707 .0707];
    end
            
    
    
    rf = [vars{c,s,t}.targets.(target_str{t})(1)*vars{c,s,t}.distances(d)/.1+ro(1),...
        vars{c,s,t}.targets.(target_str{t})(2)*vars{c,s,t}.distances(d)/.1+ro(2)];
    
    [Data{c,s,t}.time,Data{c,s,t}.x,Data{c,s,t}.y,~,...
        ~,~,~] =...
        Gen_mvt_gb(Resamp,c,s,t,ro,rf,vars{c,s,t}.time_inc);
    
%     [Data{c,s,t}.time,Data{c,s,t}.x,Data{c,s,t}.y,Data{c,s,t}.vx,...
%         Data{c,s,t}.vy,Data{c,s,t}.ax,Data{c,s,t}.ay] =...
%         Gen_mvt(Resamp,c,s,t,ro,rf,vars{c,s,t}.time_inc);
%     [Data{c,s,t}.time,Data{c,s,t}.x,Data{c,s,t}.y,Data{c,s,t}.vx,...
%         Data{c,s,t}.vy,Data{c,s,t}.ax,Data{c,s,t}.ay] ...
%         =minjerk(ro,rf,vars{c,s,t}.speeds(s),vars{c,s,t}.time_inc);
    
    Data{c,s,t}.targetposition(1)=vars{c,s,t}.targets.(target_str{t})(1)*vars{c,s,t}.distances(d)/.1+ro(1);
    Data{c,s,t}.targetposition(2)=vars{c,s,t}.targets.(target_str{t})(2)*vars{c,s,t}.distances(d)/.1+ro(2);
    
    Data{c,s,t}.startposition = ro;
    
    vars{c,s,t}.masses;
            end
        end
    end
end

parfor i=1:4*6*8
    [c,s,t] = ind2sub([4,6,8],i);
    [ shoulder1{i},elbow1{i},theta1{i},muscles1{i},act1{i},u1{i},est1{i},tnew1{i}, energy1{i}, eff_mass1{i} ] =...
        looper(Data{c,s,t}, forearm{c,s,t}, upperarm{c,s,t},vars{c,s,t},time_step,c,s,t);
    if mod(i,10)==0
        i
    end
end

for i =1:4*6*8
   
    [c,s,t] = ind2sub([4,6,4],i);
    
    shoulder{c,s,t} = shoulder1{i};
    elbow{c,s,t} = elbow1{i};
    theta{c,s,t} = theta1{i};
    eff_mass{c,s,t} = eff_mass1{i};
    muscles{c,s,t} = muscles1{i};
    act{c,s,t} = act1{i};
    u{c,s,t} = u1{i};
    est{c,s,t} = est1{i};
    tnew{c,s,t} = tnew1{i};
    energy{c,s,t} = energy1{i};
    
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

filename = sprintf('aa_%scost_03-22-2018',  minparams{L});
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
