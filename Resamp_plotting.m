
1;
%%
% clear
main_fold = 'D:\Users\Gary\Google Drive\Muscle modeling\Metabolics\'
% main_fold = 'E:\Documents\Google Drive\Muscle modeling\Metabolics\';
cd(main_fold);
if ~exist('Resamp')
    cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk_files\Data')
    load('Resamp_data_PLOTTING.mat');
    n_samp = 100;
    masses = {'0 lbs', '5 lbs', '10 lbs', '20 lbs'};
%     cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk_files\2019 Data');
%     load('aa_act2cost_09-17-2019.mat');  
%     cd(cur_fold);
end


figure(1);clf(1);
% [ha,pos] = tight_subplot(2,4,[.1 .03],[.05 .05],[.05 .05]);
Colors = parula(7);
for c = 1:4
    if c == 1 || c ==2
        speds = 1:6;
    else
        speds = 2:7;
    end
    for s = speds
        subplot(2,4,c);
        hold on
        x=.015+reshape(Resamp.T(c,s,:),1,n_samp);
        y=reshape(Resamp.P(c,s,:),1,n_samp);
        e=reshape(Resamp.P_e(c,s,:),1,n_samp);
        lo = y - e;
        hi = y + e;
        hi(hi>=.11) = 0.11;
        lo(lo<=0) = 0;
        plot(x,y,'Color',Colors(s,:));
        plot(x,hi,'Color',Colors(s,:));
        plot(x,lo,'Color',Colors(s,:));
        ylim([0 .11]);xlim([0 1.4]);
        xlabel('Time (s)');
        ylabel('Position (m)');
        subplot(2,4,c+4);
        hold on
        clear diffV
        diffV(1) = 0;
        for k = 2:100
            diffV(k) = (Resamp.P(c,s,k)-Resamp.P(c,s,k-1))/(Resamp.T(c,s,k)-Resamp.T(c,s,k-1));
        end
        diffV = sgolayfilt(diffV,1,21);
        diffV(1) = 0;
        diffV(2) = 0;
        diffV = spline([Resamp.T(c,s,1);Resamp.T(c,s,2);reshape(Resamp.T(c,s,20:20:80),4,1);Resamp.T(c,s,end-1);Resamp.T(c,s,end)],...
            [0,0,diffV(20:20:80),0,0],[Resamp.T(c,s,:)]);
        diffV = reshape(diffV,100,1);
        
        
        x=.015+reshape(Resamp.T(c,s,:),1,n_samp);
        y=(diffV)';
        e=reshape(Resamp.V_e(c,s,:),1,n_samp);
        lo = y - e;
        hi = y + e;
        hi(hi>=.7) = 0.7;
        lo(lo<0) = 0;
        plot(x,y,'Color',Colors(s,:));
        plot(x,hi,'Color',Colors(s,:));
        plot(x,lo,'Color',Colors(s,:));
        
        xlim([0 1.4]);
        ylim([0 0.7]);
        title(sprintf('%s Added Mass',masses{c}));
        xlabel('Time (s)');
        ylabel('Velocity (m)');
    end
end
speeds = {'VVVF','VVF','VF','F','M','S','VS'};

cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk_files\2019 Data');
load('aa_act2cost_09-17-2019.mat');

clear P V
for s = 1:6
    for c = 1:4
        T1 = Data{c,s,1}.time;
        P = sqrt((Data{c,s,t}.x-Data{c,s,t}.x(1)).^2+...
                (Data{c,s,t}.y-Data{c,s,t}.y(1)).^2);

        t=1;
        addpath('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files');
        muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};
        for ii = 1:length(act{c,s,t}.an)
            norm_force = vars{c,s,t}.norm_force;
            for k = 1:length(muscle_nums)
                norm_length(k) = muscles{c,s,t}.(muscle_nums{k}).length(ii)/muscles{c,s,t}.(muscle_nums{k}).l0;
                vel(k) = muscles{c,s,t}.(muscle_nums{k}).v(ii);
        %         n_f(k) = muscles{c,s,t}.(muscle_nums{k}).force(ii)/(norm_force*muscles{c,s,t}.(muscle_nums{k}).pcsa);

        %         a = est{c,s,t}.(muscle_nums{k})(ii);
                a = act{c,s,t}.(muscle_nums{k})(ii);

                check.stress(k,ii) = Fl_Fv_for(norm_length(k),vel(k),a)*norm_force;
                check.force(k,ii) = check.stress(k,ii)*muscles{c,s,t}.(muscle_nums{k}).pcsa;
            end



            A1 = [muscles{c,s,t}.an.m_arm_e(ii),...
                muscles{c,s,t}.bs.m_arm_e(ii),...
                muscles{c,s,t}.br.m_arm_e(ii),...
                0,...
                0,...
                0,...
                muscles{c,s,t}.bb.m_arm_e(ii),...
                muscles{c,s,t}.tb.m_arm_e(ii)];
            A2 = [0,...
                0,...
                0,...
                muscles{c,s,t}.da.m_arm_s(ii),...
                muscles{c,s,t}.dp.m_arm_s(ii),...
                muscles{c,s,t}.pc.m_arm_s(ii),...
                muscles{c,s,t}.bb.m_arm_s(ii),...
                muscles{c,s,t}.tb.m_arm_s(ii)];
            A=[A1;A2];

            x = [check.force(1,ii),check.force(2,ii),check.force(3,ii),...
                check.force(4,ii),check.force(5,ii),check.force(6,ii),...
                check.force(7,ii),check.force(8,ii)]';

            check.elbow{c,s,t}.torque_c(ii) = A1*x;
            check.shoulder{c,s,t}.torque_c(ii) = A2*x;

            prm.m1 = upperarm{c,s,t}.mass;
            prm.r1 = upperarm{c,s,t}.centl;
            prm.l1 = upperarm{c,s,t}.length;
            prm.i1 = upperarm{c,s,t}.Ic;

            prm.m2 = forearm{c,s,t}.mass;
            prm.r2 = forearm{c,s,t}.l_com;
            prm.r22 = forearm{c,s,t}.centl;
            prm.l2 = forearm{c,s,t}.length;
            prm.i2 = forearm{c,s,t}.Ic;

            prm.m = vars{c,s,t}.masses;

             M11 = prm.m1*prm.r1^2 + prm.i1 +...
                (prm.m+prm.m2)*(prm.l1^2+prm.r22^2+...
                2*prm.l1*prm.r22*cos(theta{c,s,t}.E(ii))) +...
                prm.i2;

            M12 = (prm.m2+prm.m)*(prm.r22^2+...
                prm.l1*prm.r22*cos(theta{c,s,t}.E(ii))) +...
                prm.i2;

            M21 = M12;

            M22 = prm.m2*prm.r2^2+prm.m*prm.l2^2+prm.i2;

            C1 = -forearm{c,s,t}.mass*forearm{c,s,t}.centl*upperarm{c,s,t}.length*(theta{c,s,t}.Ed(ii).^2).*sin(theta{c,s,t}.E(ii))-...
                2*forearm{c,s,t}.mass*forearm{c,s,t}.centl*upperarm{c,s,t}.length*theta{c,s,t}.Sd(ii).*theta{c,s,t}.Ed(ii).*sin(theta{c,s,t}.E(ii));

            C2 = forearm{c,s,t}.mass*forearm{c,s,t}.centl*upperarm{c,s,t}.length*(theta{c,s,t}.Sd(ii).^2).*sin(theta{c,s,t}.E(ii));

            qdd = [M11,M12;M21,M22] \([check.shoulder{c,s,t}.torque_c(ii); check.elbow{c,s,t}.torque_c(ii)] - [C1;C2]);

        %     qdd(1) = I \(check.shoulder{c,s,t}.torque_c(ii) - C1);
        %     qdd(1) = I \(check.elbow{c,s,t}.torque_c(ii) - C2);

            check.theta{c,s,t}.Sdd(ii) = qdd(1);
            check.theta{c,s,t}.Edd(ii) = qdd(2);
        end
        vars{c,s,t}.masses
        %% Set Init Conditions
        check.theta{c,s,t}.S(1,:) = theta{c,s,t}.S(1,:);
        check.theta{c,s,t}.E(1,:) = theta{c,s,t}.E(1,:);
        check.theta{c,s,t}.Sd(1,1:length(theta{c,s,t}.Sd(1,:))) = theta{c,s,t}.Sd(1,:);
        check.theta{c,s,t}.Ed(1,1:length(theta{c,s,t}.Ed(1,:))) = theta{c,s,t}.Ed(1,:);

        %% Calulate Kinematics
        check.theta{c,s,t}.Sd = (cumsum(check.theta{c,s,t}.Sdd)*vars{c,s,t}.time_inc+check.theta{c,s,t}.Sd(1))';
        check.theta{c,s,t}.Ed = (cumsum(check.theta{c,s,t}.Edd)*vars{c,s,t}.time_inc+check.theta{c,s,t}.Ed(1))';

        check.theta{c,s,t}.S = (cumsum(check.theta{c,s,t}.Sd)*vars{c,s,t}.time_inc+check.theta{c,s,t}.S(1));
        check.theta{c,s,t}.E = (cumsum(check.theta{c,s,t}.Ed)*vars{c,s,t}.time_inc+check.theta{c,s,t}.E(1));

        check.theta{c,s,t}.Sdd = check.theta{c,s,t}.Sdd';
        check.theta{c,s,t}.Edd = check.theta{c,s,t}.Edd';

        check.theta{c,s,t}.S = reshape(check.theta{c,s,t}.S,[length(check.theta{c,s,t}.S),1]);
        check.theta{c,s,t}.E = reshape(check.theta{c,s,t}.E,[length(check.theta{c,s,t}.E),1]);

        check.theta{c,s,t}.Sd = reshape(check.theta{c,s,t}.Sd,[length(check.theta{c,s,t}.Sd),1]);
        check.theta{c,s,t}.Ed = reshape(check.theta{c,s,t}.Ed,[length(check.theta{c,s,t}.Ed),1]);

        check.theta{c,s,t}.Sdd = reshape(check.theta{c,s,t}.Sdd,[length(check.theta{c,s,t}.Sdd),1]);
        check.theta{c,s,t}.Edd = reshape(check.theta{c,s,t}.Edd,[length(check.theta{c,s,t}.Edd),1]);


        %% Calc Error

        s_error(c,s,t,d) = sum(abs((check.theta{c,s,t}.S - theta{c,s,t}.S)));
        e_error(c,s,t,d) = sum(abs((check.theta{c,s,t}.E - theta{c,s,t}.E)));

        x = cos(check.theta{c,s,t}.S)*upperarm{c,s,t}.length +...
            forearm{c,s,t}.length*cos(check.theta{c,s,t}.S+check.theta{c,s,t}.E);
        y = sin(check.theta{c,s,t}.S)*upperarm{c,s,t}.length +...
            forearm{c,s,t}.length*sin(check.theta{c,s,t}.S+check.theta{c,s,t}.E);

        P = sqrt((x-Data{c,s,t}.x(1)).^2+(y-Data{c,s,t}.y(1)).^2);

        V = [0;diff(P)]/0.0025;
        if c == 4 || c == 3
            pltcolor = Colors(s+1,:);
        else
            pltcolor = Colors(s,:);
        end
        subplot(2,4,c);
        %         axes(ha(c));
        a=gca;
        set(a,'tickdir','out');
        set(a,'XtickLabel',a.XTick);
        set(a,'Ytick',[0,.02,.04,.06,.08,.10],'YtickLabel',[0,2,4,6,8,10]);
        %         subplot(2,4,c);
        hold on
        plot(T1,P,'--','Color',pltcolor,'linewidth',2);
        ylim([0 .11]);xlim([0 1.4]);
        title(sprintf('%s Added Mass',masses{c}));
        xlabel('Time (s)');
        ylabel('Position (cm)');
        subplot(2,4,c+4);
        %         axes(ha(c+4));
        a=gca;
        set(a,'tickdir','out');
        set(a,'XtickLabel',a.XTick);
        set(a,'YtickLabel',a.YTick);
        %         subplot(2,4,c+4);
        hold on
        plot(T1,V,'--','Color',pltcolor,'linewidth',2);
        xlim([0 1.4]);
        title(sprintf('%s Added Mass',masses{c}));
        xlabel('Time (s)');
        ylabel('Velocity (m/s)');
    end
end
% subplot(2,4,1);
% legend(plots,{'VVVF','VVF','VF','F','M','S','VS'});
% beautifyfig;

diffV(1) = 0
for k = 2:100
    diffV(k) = (Resamp.P(1,1,k)-Resamp.P(1,1,k-1))/(Resamp.T(1,1,k)-Resamp.T(1,1,k-1));
end
