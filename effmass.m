time_step = .005;
tic

target_str = {'one','two','thr','fou'};
input_normforce = 31.84E4;
% masses = [0,5,10,20]/2.2;

clearvars -except time_step target_str input_normforce masses minparams L
for c=1:4
    for s=1:6
        for t = 1:4
            for d = 2:2
                
    vars{c,s,t,d}.time_inc = time_step;
%     vars{c,s,t,d}.masses = [0,3,5,8]/2.2;
    vars{c,s,t,d}.masses = [0,5,10,20]/2.2;
%     vars{c,s,t,d}.masses = masses/2.2;
    % vars{c,s,t,d}.speeds = {'VS','S','M','F','VF','VVF'};
    vars{c,s,t,d}.speeds = [.45, .55, .70, .85, 1.00, 1.15];
    switch c
        case 1
        vars{c,s,t,d}.speeds = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
        case 2
        vars{c,s,t,d}.speeds = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
        case 3
        vars{c,s,t,d}.speeds = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
        case 4
        vars{c,s,t,d}.speeds = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
    end
    vars{c,s,t,d}.distances = [.05,.1,.15,.2];
    vars{c,s,t,d}.targets.one = [.0707 .0707];
    vars{c,s,t,d}.targets.two = [-.0707 .0707];
    vars{c,s,t,d}.targets.thr = [-.0707 -.0707];
    vars{c,s,t,d}.targets.fou = [.0707 -.0707];
    vars{c,s,t,d}.norm_force = input_normforce;%200E4;%31.8E4;
%     vars{c,s,t,d}.minparam = minparams{L};
    
    upperarm{c,s,t,d}.length = .33; % meters
    upperarm{c,s,t,d}.l_com = upperarm{c,s,t,d}.length/2;
    upperarm{c,s,t,d}.centl = upperarm{c,s,t,d}.length/2;
    upperarm{c,s,t,d}.mass = 1.93284; %kg
%     upperarm{c,s,t,d}.Ic = (1/12)*upperarm{c,s,t,d}.mass*upperarm{c,s,t,d}.length.^2;
    upperarm{c,s,t,d}.Ic = .0141;
    
    forearm{c,s,t,d}.mass = 1.5186;%;+2*vars{c,s,t,d}.masses(c);
    forearm{c,s,t,d}.length = upperarm{c,s,t,d}.length*1.3; % meters, taken from An iterative optimal control and estimation design for nonlinear stochastic system
    forearm{c,s,t,d}.l_com = forearm{c,s,t,d}.length *2/3;
    [forearm{c,s,t,d}] = calc_forearmI(forearm{c,s,t,d},vars{c,s,t,d}.masses(c));
    forearm{c,s,t,d};
%     forearm{c,s,t,d}.Ic = .01882;


    shoulder{c,s,t,d} = [];
    elbow{c,s,t,d} = [];
    theta{c,s,t,d} = [];
    
%     ro = [0,.4];
    ro = [-.0758,0.4878];
    rf = [vars{c,s,t,d}.targets.(target_str{t})(1)*vars{c,s,t,d}.distances(d)/.1+ro(1),...
        vars{c,s,t,d}.targets.(target_str{t})(2)*vars{c,s,t,d}.distances(d)/.1+ro(2)];
    
    [Data{c,s,t,d}.time,Data{c,s,t,d}.x,Data{c,s,t,d}.y,Data{c,s,t,d}.vx,...
        Data{c,s,t,d}.vy,Data{c,s,t,d}.ax,Data{c,s,t,d}.ayy] ...
        =minjerk(ro,rf,vars{c,s,t,d}.speeds(s),vars{c,s,t,d}.time_inc);
    
    Data{c,s,t,d}.targetposition(1)=vars{c,s,t,d}.targets.(target_str{t})(1)*vars{c,s,t,d}.distances(d)/.1+ro(1);
    Data{c,s,t,d}.targetposition(2)=vars{c,s,t,d}.targets.(target_str{t})(2)*vars{c,s,t,d}.distances(d)/.1+ro(2);
    
    Data{c,s,t,d}.startposition = ro;
    
    vars{c,s,t,d}.masses;
    
    [ shoulder{c,s,t,d} , elbow{c,s,t,d} , theta{c,s,t,d} , eff_mass{c,s,t,d}] = ...
        Calc_kine_test( Data{c,s,t,d}, forearm{c,s,t,d}, upperarm{c,s,t,d},vars{c,s,t,d}.time_inc,vars{c,s,t,d}.masses(c));
            end
        end
    end
end