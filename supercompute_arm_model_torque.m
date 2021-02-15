clc
time_step = .005;

target_str = {'one','two','thr','fou','fiv','six','sev','eig'};
input_normforce = 200E4;
minparams = {'umberger','stress','force','act','drive','stress2','force2','act2','drive2'};

t_act= 0.05;%.0050;
t_deact = 0.066;%.066;

fold = pwd;
% load('Data\Resamp_data_gb2.mat');
% 
if ~exist('L')
    L = 2;
end
if ~exist('num_work')
    num_work = 1;
end

if num_work > 1
    mycluster=parcluster('local');
    mycluster.NumWorkers=num_work;
    parpool(num_work);
end
% L=5
% for L= [1,2,3,5,6,7,8,9]
tic
minfail = [];
if exist('rnjesus');
    if rnjesus
        filename = sprintf('aa_%scost_7-17-2020_rng',  minparams{L});
    else
        filename = sprintf('aa_%scost_7-17-2020',  minparams{L});
    end
else
    rnjesus = 0;
    filename = sprintf('aa_%scost_7-17-2020',  minparams{L});
end

fprintf('Processing %s \n',filename);
clearvars -except time_step target_str input_normforce masses minparams L Resamp tic fold minfail rnjesus

rng(57);
speeds = 1;
masses = 1;
% subj_mass = rand(10,1)*60+130;
% subj_heights = rand(10,1)*30+150;
subj_mass = 60*(2.2);
subj_heights = 175/100;

upperarm_mass = 0.028*subj_mass/2.2;
lowerarm_mass = 0.022*subj_mass/2.2;
upperarm_length = (.825-.632)*subj_heights;
lowerarm_length = (.632-0.425)*subj_heights; %hand is half the length
upperarm_lcom = 0.436*upperarm_length;
lowerarm_lcom = 0.682*lowerarm_length;

for c=1:length(masses)
    for s=1:length(speeds)
        for t = 1:8
            d=2;
% if sum(isnan(Resamp.P(c,s,t,:)))>10 || sum(Resamp.P(c,s,t,:))<.1
%     continue
% end
vars{c,s,t}.rnjesus = rnjesus;
vars{c,s,t}.time_inc = time_step;
vars{c,s,t}.speeds = speeds;% squeeze(Resamp.Endt(c,:,t));
vars{c,s,t}.masses = masses(c);

vars{c,s,t}.distances = [.05,.1,.15,.2];
vars{c,s,t}.norm_force = input_normforce;%200E4;%31.8E4;
vars{c,s,t}.minparam = minparams{L};
vars{c,s,t}.L = L;

% http://www.kdm.p.lodz.pl/articles/2017/3/21_3_4.pdf
% https://www.ele.uri.edu/faculty/vetter/BME207/anthropometric-data.pdf
forearm{c,s,t}.mass = lowerarm_mass;%;+2*vars{c,s,t}.masses(c);
forearm{c,s,t}.length = lowerarm_length; % meters, taken from An iterative optimal control and estimation design for nonlinear stochastic system
forearm{c,s,t}.l_com = 0.682*lowerarm_length;

forearm{c,s,t}.l_com = (.417*(.632-.480)*subj_heights*.016*subj_mass +...
        ((.632-.480)*subj_heights+.515*((.480-.370)/2)*subj_heights)*.006*subj_mass)/...
        (.016*subj_mass+.006*subj_mass); % Using Enoka
    
[forearm{c,s,t}] = calc_forearmI(forearm{c,s,t},vars{c,s,t}.masses);

upperarm{c,s,t}.length = upperarm_length; % meters
upperarm{c,s,t}.l_com = 0.436*upperarm_length;
upperarm{c,s,t}.centl = 0.436*upperarm_length;
upperarm{c,s,t}.mass = upperarm_mass; %kg
upperarm{c,s,t}.Ic = .0141;

shoulder{c,s,t} = [];
elbow{c,s,t} = [];
theta{c,s,t} = [];

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

% % This is for generated movement from resamp
% [Data{c,s,t}] = Gen_mvt_gb(Resamp,c,s,t,ro,rf,vars{c,s,t}.time_inc);
% 
% Data{c,s,t}.targetposition(1)=vars{c,s,t}.targets.(target_str{t})(1)*vars{c,s,t}.distances(d)/.1+ro(1);
% Data{c,s,t}.targetposition(2)=vars{c,s,t}.targets.(target_str{t})(2)*vars{c,s,t}.distances(d)/.1+ro(2);
% Data{c,s,t}.startposition = ro;

% This is for minjerk
Data{c,s,t} = minjerk(ro,rf,vars{c,s,t}.speeds(s),vars{c,s,t}.time_inc);
Data{c,s,t}.targetposition = rf;
Data{c,s,t}.startposition  = ro;
Data{c,s,t}.added_mass = masses(c);
Data{c,s,t}.move_dur = speeds(s);

vars{c,s,t}.masses;
        end
    end
end

parfor_progress(length(masses)*length(speeds)*8);
for i=1:length(masses)*length(speeds)*8
    [c,s,t] = ind2sub([length(masses),length(speeds),8],i);
    if isempty(Data{c,s,t})
        continue
    end
    [ shoulder1{i},...
        elbow1{i},...
        theta1{i},...
        eff_mass1{i}] =...
            looper2_torque(Data{c,s,t},...
            forearm{c,s,t},...
            upperarm{c,s,t},...
            vars{c,s,t});
   parfor_progress;
end
parfor_progress(0);

for i = 1:length(masses)*length(speeds)*8
    [c,s,t] = ind2sub([length(masses),length(speeds),8],i);
    
    shoulder{c,s,t} = shoulder1{i};
    elbow{c,s,t} = elbow1{i};
    theta{c,s,t} = theta1{i};
    eff_mass{c,s,t} = eff_mass1{i};
    Data{c,s,t}.theta = theta1{i};
    Data{c,s,t}.shoulder = shoulder1{i};
    Data{c,s,t}.elbow = elbow1{i};
    Data{c,s,t}.eff_mass = eff_mass1{i};
    Data{c,s,t}.torque2 = time_step*(sum(Data{c,s,t}.elbow.torque.^2)+sum(Data{c,s,t}.shoulder.torque.^2));
end
% Checkers_test;
clear shoulder1 elbow1 theta1 muscles1 act1 est1 tnew1
toc

% filename = sprintf('torque_7-17',  minparams{L});
% save(filename);
% if exist(strcat(filename,'.mat'),'file')==2
%     fprintf('Complete. Saved as: %s\n',filename');
% end 

%% Write the data
fprintf('Writing data to excel.\n')
excel_file = 'sum_torque2.csv';
fileID=fopen(excel_file,'w');

delcount=0;

A={'movedur,'...
   'added_mass,'...
   'target,'...
   'effmass,'...
   'effmass_mean,'...
   'sum_t2'...
   };

fprintf(fileID,'%s',A{:});
fprintf(fileID,'\n');

for c=1:length(masses)
    for s=1:length(speeds)
       for t = 1:8
           data = Data{c,s,t};
            if ~isempty(data)
                A={data.move_dur,...
                   data.added_mass,...
                   t,...
                   data.eff_mass(1),...
                   mean(data.eff_mass),...
                   data.torque2};
                   dlmwrite(excel_file,A,'-append')
            end
        end
    end
end

fclose(fileID);

% end
