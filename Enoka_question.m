        
 %% Init vars
% 
% cd('D:\Users\Gary\OneDrive\Muscle modeling\Metabolics Files');
% if exist('Data')~=1
%     load('XY_outstop_data.mat');
%     fprintf('Loaded data \n');
% end

%% Init Vars
clear
% global time_step
time_step = .005;
tic

target_str = {'one','two','thr','fou'};

clear x forearm upperarm muscles A1 A2 Aeq beq
for c=1:4
%     for s=1:6
        for t = 1:4
%             for d = 1:4
                
    vars{c,t}.time_inc = time_step;
    vars{c,t}.masses = [0,3,5,8]/2.2;
    % vars{c,t}.speeds = {'VS','S','M','F','VF','VVF'};
    vars{c,t}.speeds = [.95,1.0098,1.0395,1.085];
    vars{c,t}.distances = [.1];
    vars{c,t}.targets.one = [.0707 .0707];
    vars{c,t}.targets.two = [-.0707 .0707];
    vars{c,t}.targets.thr = [-.0707 -.0707];
    vars{c,t}.targets.fou = [.0707 -.0707];
    
    forearm{c,t}.mass = 1+2*vars{c,t}.masses(c);
    forearm{c,t}.length = .30; % meters, taken from An iterative optimal control and estimation design for nonlinear stochastic system
    forearm{c,t}.centl = .11;
    forearm{c,t}.Ic = (1/3)*forearm{c,t}.mass*forearm{c,t}.length^2;

    upperarm{c,t}.length = .33; % meters
    upperarm{c,t}.centl = .16;
    upperarm{c,t}.mass = 1.4; %kg
    upperarm{c,t}.Ic = (1/3)*upperarm{c,t}.mass*upperarm{c,t}.length^2;

    shoulder{c,t} = [];
    elbow{c,t} = [];
    theta{c,t} = [];

    [Data{c,t}.time,Data{c,t}.x,Data{c,t}.y,Data{c,t}.vx,...
        Data{c,t}.vy,Data{c,t}.ax,Data{c,t}.ayy] ...
        =minjerk([0,.4],[vars{c,t}.targets.(target_str{t})(1)*vars{c,t}.distances/.1,...
        vars{c,t}.targets.(target_str{t})(2)*vars{c,t}.distances/.1+.4],vars{c,t}.speeds(c),vars{c,t}.time_inc);
    
    Data{c,t}.targetposition(1)=vars{c,t}.targets.(target_str{t})(1)*vars{c,t}.distances/.1;
    Data{c,t}.targetposition(2)=vars{c,t}.targets.(target_str{t})(2)*vars{c,t}.distances/.1+.4;
    
    Data{c,t}.startposition = [0,.4];
%             end
        end
%     end
end

parfor i=1:4*4
    [c,t] = ind2sub([4,4],i);
    [ shoulder1{i},elbow1{i},theta1{i},muscles1{i},act1{i},u1{i},est1{i},tnew1{i} ] =...
        looper(Data{c,t}, forearm{c,t}, upperarm{c,t},vars{c,t},time_step);
i  
end
% parfor_progress(0);
toc

for i = 1:4*4
   
    [c,t] = ind2sub([4,4],i);
    
    shoulder{c,t} = shoulder1{i};
    elbow{c,t} = elbow1{i};
    theta{c,t} = theta1{i};
    muscles{c,t} = muscles1{i};
    act{c,t} = act1{i};
    u{c,t} = u1{i};
    est{c,t} = est1{i};
    tnew{c,t} = tnew1{i};
    
    shoulder1{i} = [];
    elbow1{i} = [];
    theta1{i} = [];
    muscles1{i} = [];
    act1{i} = [];
    est1{i} = [];
    tnew1{i} = [];
    
end
clear shoulder1 elbow1 theta1 muscles1 act1 est1 tnew1

filename = 'test_72017_140force_cap';
% save(filename);
if exist(strcat(filename,'.mat'),'file')==2
    fprintf('File saved as %s \n',filename');
end

%% Plottign
muscle_nums = {'one','two','thr','fou','fiv','six'};
ColorSet = parula(4);

% mass_legend = {'0 lbs','3 lbs','5 lbs','8 lbs'};
mass_legend = {'0 lbs','3 lbs','5 lbs','8 lbs'};

%% Neural Activation
figure(1);clf(1); hold on;
subplot(2,2,1);hold on;
for c = 1:4
        
        usum=0;
        for t=1:4
            for k = 1:length(muscle_nums)
                usum =  usum + sum(abs(u{c,t}.(muscle_nums{k})))*time_step;
                if ~isreal(usum)
                    1;
                end
            end
        end
        usum1(c) = usum;
        ColorPref(c) = plot(vars{c,t}.speeds(c),usum,'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        

end
plot(vars{c,t}.speeds,usum1,'k-');
legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel('Neural Drive');

%% Muscle Active State
figure(1);%clf(1); hold on;
subplot(2,2,2);hold on;
for c = 1:4
        
    usum=0;
    for t=1:4
        for k = 1:length(muscle_nums)
            ind = act{c,t}.(muscle_nums{k})>1E-6;
            if sum(ind)>0
                usum =  usum + mean(act{c,t}.(muscle_nums{k}))*time_step;
            end
            if ~isreal(usum)
                1;
            end
        end
    end
    usum1(c) = usum;
    ColorPref(c) = plot(vars{c,t}.speeds(c),usum,'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));

end

    plot(vars{c,t}.speeds,usum1,'k-');
legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel('Active State');


%% Total Force
%figure(2);clf(2); hold on;
subplot(2,2,3);hold on;
for c = 1:4
        
        fsum=0;
        for t=1:4
            for k = 1:length(muscle_nums)
                fsum =  fsum + sum(muscles{c,t}.(muscle_nums{k}).force)*time_step;
            end
        end
        fsum1(c) = fsum;
        ColorPref(c) = plot(vars{c,t}.speeds(c),fsum,'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        
end

    plot(vars{c,t}.speeds,fsum1,'k-');
legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel('Muscle Force');

%% Joint Torque
%figure(4);clf(4); hold on;
subplot(2,2,4);hold on;
for c = 1:4
        
        tsum=0;
        for t=1:4
            tsum =  tsum + (sum(elbow{c,t}.torque)+sum(shoulder{c,t}.torque))*time_step;
        end
        tsum1(c) = tsum;
        ColorPref(c) = plot(vars{c,t}.speeds(c),tsum,'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        
end

    plot(vars{c,t}.speeds,tsum1,'k-');
legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel('Joint Torque');

cd('d:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files\Graphs');
% savefig('Averages_normforce');
cd('d:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files');


