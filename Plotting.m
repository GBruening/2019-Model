
muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};
ColorSet = parula(5);
ColorSet = ColorSet(1:4,:);

% mass_legend = {'0 lbs','3 lbs','5 lbs','8 lbs'};
mass_legend = {'0 lbs','5 lbs','10 lbs','20 lbs'};

d=2;
if ~exist('expo')
    expo = 2;
end

n_tar = 8;

%% Neural Activation
figure(1);clf(1); hold on;
subplot(2,2,1);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            usum=0;
            for t=1:n_tar
                for k = 1:length(muscle_nums)
                    sum(abs(u{c,subj,s,t}.(muscle_nums{k}).^expo))*(time_step/10);
                    usum =  usum + sum(abs(u{c,subj,s,t}.(muscle_nums{k}).^expo))*(time_step)*muscles{c,subj,s,t}.(muscle_nums{k}).m;
    %                 fprintf('YOU DIDN"T FIX THE DRIVE TIME STEP INTEGRAL DUMBO');
                    if ~isreal(usum)
                        1;
                    end
                end
            end
            usum1(c,subj,s) = usum/vars{c,subj,s,t}.speeds(s)/n_tar;
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),usum1(c,subj,s)/vars{c,subj,s,t}.speeds(s)/n_tar,...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
end

legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel(sprintf('Neural Drive Dot^%g',expo));

drawnow;pause(.1);

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(1);savefig('Neural Dirve');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');%% Neural Activation

%% Average Neural Activation
% figure(1);%clf(1); hold on;
% subplot(2,2,2);hold on;
% for c = 1:4
%     for s = 1:6
%         
%         usum=0;
%         for t=1:4
%             for k = 1:length(muscle_nums)
%                 ind = u{c,subj,s,t}.(muscle_nums{k})~=0;
%                 if sum(ind)>0
%                     usum =  usum + mean(abs(u{c,subj,s,t}.(muscle_nums{k})(ind)))*time_step;
%                 end
%                 if ~isreal(usum)
%                     1;
%                 end
%             end
%         end
%         usum1(c,subj,s) = usum/vars{c,subj,s,t}.speeds(s);
%         ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),usum/vars{c,subj,s,t}.speeds(s),'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
%         
%     end
%     plot(vars{c,subj,s,t}.speeds,usum1(c,:),'k-');
% end
% 
% legend([ColorPref(:)],mass_legend);
% xlabel('Movement Duration');ylabel('Neural Drive Avg Dot');

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(1);savefig('Neural Dirve');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');


%% Muscle Active State
figure(1);%clf(1); hold on;
subplot(2,2,2);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            asum=0;
            for t=1:n_tar
                for k = 1:length(muscle_nums)
                    ind = act{c,subj,s,t}.(muscle_nums{k})>1E-9;
                    if sum(ind)>0
                        asum =  asum + sum(act{c,subj,s,t}.(muscle_nums{k})(ind).^expo)*time_step;
                    end
                    if ~isreal(asum)
                        1;
                    end
                end
            end
            asum1(c,subj,s) = asum/vars{c,subj,s,t}.speeds(s)/n_tar;
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),asum/vars{c,subj,s,t}.speeds(s)/n_tar,...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,asum1(c,:),'k-');
end

legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel(sprintf('Active State Dot^%g',expo));


%% Total Muscle Force
%figure(2);clf(2); hold on;
subplot(2,2,3);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            fsum=0;
            for t=1:n_tar
                for k = 1:length(muscle_nums)
                    fsum =  fsum + sum(muscles{c,subj,s,t}.(muscle_nums{k}).force.^expo)*time_step;
                end
            end
            fsum1(c,subj,s) = fsum/(vars{c,subj,s,t}.speeds(s)/n_tar);
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),fsum/vars{c,subj,s,t}.speeds(s)/n_tar,...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,fsum1(c,:),'k-');
end

legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel(sprintf('Muscle Force Dot^%g',expo));

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(2);savefig('Total Force');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');

%% Total Output Force
figure(2);clf(2); hold on;
% subplot(2,2,3);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            f_out = 0;
            for t=1:n_tar
                fx_1 = -sin(theta{c,subj,s,t}.S).*...
                    shoulder{c,subj,s,t}.torque/...
                    upperarm{c,subj,s,t}.length;
                fy_1 = cos(theta{c,subj,s,t}.S).*...
                    shoulder{c,subj,s,t}.torque/...
                    upperarm{c,subj,s,t}.length;
                fx_2 = -sin(theta{c,subj,s,t}.S+theta{c,subj,s,t}.E).*...
                    elbow{c,subj,s,t}.torque/...
                    forearm{c,subj,s,t}.length + fx_1;
                fy_2 = cos(theta{c,subj,s,t}.S+theta{c,subj,s,t}.E).*...
                    elbow{c,subj,s,t}.torque/...
                    forearm{c,subj,s,t}.length + fy_1;
                fx_2 = sum(fx_2)*0.0025;
                fy_2 = sum(fy_2)*0.0025;
                f_out = f_out + sqrt(fx_2.^2+fy_2.^2);
            end
            f_out1(c,subj,s) = f_out/(vars{c,subj,s,t}.speeds(s)/n_tar);
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),f_out/vars{c,subj,s,t}.speeds(s)/n_tar,...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,fsum1(c,:),'k-');
end
1;
% legend([ColorPref(:)],mass_legend);
% xlabel('Movement Duration');ylabel(sprintf('Muscle Force Dot^%g',expo));

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(2);savefig('Total Force');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');

%% Muscle Stress
%figure(3);clf(3); hold on;
% subplot(2,2,3);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end     
            ssum=0; 
            for t=1:4
                for k = 1:length(muscle_nums)
                    ssum =  ssum + sum(muscles{c,subj,s,t}.(muscle_nums{k}).force/muscles{c,subj,s,t}.(muscle_nums{k}).pcsa)*time_step;
                end
            end
            ssum1(c,subj,s) = ssum/vars{c,subj,s,t}.speeds(s)/n_tar;
%             ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),ssum/vars{c,subj,s,t}.speeds(s),'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,ssum1(c,:),'k-');
end

% legend([ColorPref(:)],mass_legend);
% xlabel('Movement Duration');ylabel('Muscle Stress Dot');
% 
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(3);savefig('Muscle Stress');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');
%% Joint Torque
%figure(4);clf(4); hold on;
subplot(2,2,4);hold on;
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
        
            tsum=0;
            for t=1:n_tar
                tsum =  tsum + (sum(abs(elbow{c,subj,s,t}.torque.^expo))...
                    +sum(abs(shoulder{c,subj,s,t}.torque.^expo)))*time_step;
            end
            tsum1(c,subj,s) = tsum/vars{c,subj,s,t}.speeds(s)/n_tar;
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),...
                tsum/vars{c,subj,s,t}.speeds(s)/n_tar,'o','Color',ColorSet(c,:),...
                'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,tsum1(c,:),'k-');
end

legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel(sprintf('Joint Torque Dot^%g',expo));

%% Umberger
figure(44);clf(44); hold on;
% subplot(2,2,3);hold on;
energies = {'h_SL','h_AM','w'};
for c = 1:4
    for subj = 1:8
        for s = 1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            bsum=0;
            for t=1:n_tar
                for k = 1:length(muscle_nums)
                    for p1 = 1:3
                        bsum =  bsum + sum(energy{c,subj,s,t}.(muscle_nums{k}).(energies{p1}).^expo);%*time_step;
                    end
                end
                bsum;
            end
            bsum1(c,subj,s) = bsum/vars{c,subj,s,t}.speeds(s)/n_tar;
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),bsum/vars{c,subj,s,t}.speeds(s)/n_tar,...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
    end
%     plot(vars{c,subj,s,t}.speeds,bsum1(c,:),'k-');
end

legend([ColorPref(:)],mass_legend);
xlabel('Movement Duration');ylabel(sprintf('Umberger Model Dot^%g',expo));
            
drawnow;
pause(5);
%% Writing Data

mp_data = readtable('met_power_data.csv');
eff_masses = [2.47,4.73,6.99,11.5];
met_data = zeros([4,8,7]);
for c=1:4
    for subj = 1:8
        if subj==7
            1;
        end
        tempdurs = mp_data.movedur(mp_data.subject==subj & mp_data.effmass==eff_masses(c));
        tempmet = mp_data.metpowernet(mp_data.subject==subj & mp_data.effmass==eff_masses(c));
        counter = 0;
        for s=1:7
            if vars{c,subj,3,t}.speeds(s)==0
                continue
            end
            counter = counter + 1;
            met_data(c,subj,s) = tempmet(counter);
            movedur(c,subj,s) = tempdurs(counter);
        end
%         for k = 1:length(tempdurs)
%             [~,s]=min(abs(tempdurs(k)-vars{c,subj,3,t}.speeds));
%             if s==1 && subj==7 && c == 1
%                 1;
%             end
%             if isempty(vars{c,subj,s,t})
%                 continue
%             end
%             met_data(c,subj,s) = tempmet(k);
%             if met_data(c,subj,s)==0
%                 1;
%             end
%         end
    end
end

for c=1:4
    for subj=1:8
        for s=1:7
            if isempty(vars{c,subj,s,t})
                continue
            end
            if met_data(c,subj,s)==0
                1;
            end
            A={c,subj,s,movedur(c,subj,s),met_data(c,subj,s),k1,...
                tsum1(c,subj,s),f_out1(c,subj,s),fsum1(c,subj,s),...
                ssum1(c,subj,s),asum1(c,subj,s),usum1(c,subj,s),...
                bsum1(c,subj,s)};
            dlmwrite(excel_file,A,'-append');
        end
    end
end
1;
%% Save files

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% figure(4);savefig('Joint Torque');
% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files');

% cd('C:\Users\garrick\SkyDrive\Muscle modeling\Min_jerk files\Graphs');
% cd('d:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files\Graphs');
% savefig('Averages_veryhighnormforce');
% cd('d:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files');

