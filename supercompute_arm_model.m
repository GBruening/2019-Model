clc
time_step = .0025;

target_str = {'one','two','thr','fou','fiv','six','sev','eig'};
input_normforce = 200E4;
minparams = {'umberger','stress','force','act','drive','stress2','force2','act2','drive2'};

t_act= 0.05;%.0050;
t_deact = 0.066;%.066;

fold = pwd;
load('Resamp_data_gb2.mat');

% mycluster=parcluster('local');
% mycluster.NumWorkers=7;
parpool(24);
% L=5
% for L=2:2
tic
minfail = [];
if exist('rnjesus');
    if rnjesus
        filename = sprintf('aa_%scost_01-21-2019_rng',  minparams{L});
    else
        filename = sprintf('aa_%scost_01-21-2019',  minparams{L});
    end
else
    rnjesus = 0;
    filename = sprintf('aa_%scost_01-21-2019_rng',  minparams{L});
end

fprintf('Processing %s \n',filename);
clearvars -except time_step target_str input_normforce masses minparams L Resamp tic fold minfail rnjesus

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
                d=2;
if sum(isnan(Resamp.P(c,subj,s,t,:)))>10 || sum(Resamp.P(c,subj,s,t,:))<.1
    continue
end
vars{c,subj,s,t}.rnjesus = rnjesus;
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
vars{c,subj,s,t}.L = L;

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

[Data{c,subj,s,t}] = Gen_mvt_gb(Resamp,c,subj,s,t,ro,rf,vars{c,subj,s,t}.time_inc);

Data{c,subj,s,t}.targetposition(1)=vars{c,subj,s,t}.targets.(target_str{t})(1)*vars{c,subj,s,t}.distances(d)/.1+ro(1);
Data{c,subj,s,t}.targetposition(2)=vars{c,subj,s,t}.targets.(target_str{t})(2)*vars{c,subj,s,t}.distances(d)/.1+ro(2);

Data{c,subj,s,t}.startposition = ro;

vars{c,subj,s,t}.masses;
            end
        end
    end
end

parfor_progress(4*8*7*8);
for i=1:4*8*7*8
    [c,subj,s,t] = ind2sub([4,8,7,8],i);
    if isempty(Data{c,subj,s,t})
        continue
    end
    [ shoulder1{i},elbow1{i},theta1{i},muscles1{i},act1{i},u1{i},est1{i},tnew1{i}, energy1{i}, eff_mass1{i} ] =...
        looper2(Data{c,subj,s,t}, forearm{c,subj,s,t}, upperarm{c,subj,s,t},vars{c,subj,s,t},time_step,c,subj,s,t);
    parfor_progress;
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
% Checkers_test;
clear shoulder1 elbow1 theta1 muscles1 act1 est1 tnew1

toc
% fprintf('aa_%scost_01-21-2019',  minparams{L});
if exist('rnjesus');
    if rnjesus
        filename = sprintf('aa_%scost_01-21-2019_rng',  minparams{L});
    else
        filename = sprintf('aa_%scost_01-21-2019',  minparams{L});
    end
else
    filename = sprintf('aa_%scost_01-21-2019_rng',  minparams{L});
end

save(filename);
if exist(strcat(filename,'.mat'),'file')==2
    fprintf('Complete. Saved as: %s\n',filename');
end 

% end
