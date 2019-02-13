muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};

loop = 1;
v1 = zeros(4*5*4,1000);
c1 = zeros(4*5*4);
s1 = zeros(4*5*4);
t1 = zeros(4*5*4);
figure(2);clf(2);


for s=2:6
    for t = 1:4
        for c=1:2
for ii = 1:length(act{1,s,t}.an)
    norm_force = vars{1,s,t}.norm_force;
    for k = 1:length(muscle_nums)
        norm_length(k) = muscles{1,s,t}.(muscle_nums{k}).length(ii)/muscles{1,s,t}.(muscle_nums{k}).l0;
        vel(k) = muscles{1,s,t}.(muscle_nums{k}).v(ii);
        switch c
            case 1
                a = act{1,s,t}.(muscle_nums{k})(ii);
            case 2
                a = act{1,s,t}.(muscle_nums{k})(ii);
            case 3
                a = act{1,s-1,t}.(muscle_nums{k})(ii);
            case 4
                a = act{1,s-1,t}.(muscle_nums{k})(ii);                
        end
        
        check.stress(k,ii) = Fl_Fv_for(norm_length(k),vel(k),a)*norm_force;
        check.force(k,ii) = check.stress(k,ii)*muscles{1,s,t}.(muscle_nums{k}).pcsa;
    end
    
    A1 = [muscles{1,s,t}.an.m_arm_e(ii),...
        muscles{1,s,t}.bs.m_arm_e(ii),...
        muscles{1,s,t}.br.m_arm_e(ii),...
        0,...
        0,...
        0,...
        muscles{1,s,t}.bb.m_arm_e(ii),...
        muscles{1,s,t}.tb.m_arm_e(ii)];
    A2 = [0,...
        0,...
        0,...
        muscles{1,s,t}.da.m_arm_s(ii),...
        muscles{1,s,t}.dp.m_arm_s(ii),...
        muscles{1,s,t}.pc.m_arm_s(ii),...
        muscles{1,s,t}.bb.m_arm_s(ii),...
        muscles{1,s,t}.tb.m_arm_s(ii)];
    A=[A1;A2];
    
    x = [check.force(1,ii),check.force(2,ii),check.force(3,ii),...
        check.force(4,ii),check.force(5,ii),check.force(6,ii),...
        check.force(7,ii),check.force(8,ii)]';

    check.elbow{1,s,t}.torque_c(ii) = A1*x;
    check.shoulder{1,s,t}.torque_c(ii) = A2*x;
    
    prm.m1 = upperarm{1,s,t}.mass;
    prm.r1 = upperarm{1,s,t}.centl;
    prm.l1 = upperarm{1,s,t}.length;
    prm.i1 = upperarm{1,s,t}.Ic;
    
    prm.m2 = forearm{1,s,t}.mass;
    prm.r2 = forearm{1,s,t}.l_com;
    prm.r22 = forearm{1,s,t}.centl;
    prm.l2 = forearm{1,s,t}.length;
    prm.i2 = forearm{1,s,t}.Ic;
    
    prm.m = vars{1,s,t}.masses;
    
    M11 = prm.m1*prm.r1^2 + prm.i1 +...
        (prm.m+prm.m2)*(prm.l1^2+prm.r22^2+...
        2*prm.l1*prm.r22*cos(theta{1,s,t}.E(ii))) +...
        prm.i2;

    M12 = (prm.m2+prm.m)*(prm.r22^2+...
        prm.l1*prm.r22*cos(theta{1,s,t}.E(ii))) +...
        prm.i2;

    M21 = M12;

    M22 = prm.m2*prm.r2^2+prm.m*prm.l2^2+prm.i2;
    
    C1 = -forearm{1,s,t}.mass*forearm{1,s,t}.centl*upperarm{1,s,t}.length*(theta{1,s,t}.Ed(ii).^2).*sin(theta{1,s,t}.E(ii))-...
        2*forearm{1,s,t}.mass*forearm{1,s,t}.centl*upperarm{1,s,t}.length*theta{1,s,t}.Sd(ii).*theta{1,s,t}.Ed(ii).*sin(theta{1,s,t}.E(ii));

    C2 = forearm{1,s,t}.mass*forearm{1,s,t}.centl*upperarm{1,s,t}.length*(theta{1,s,t}.Sd(ii).^2).*sin(theta{1,s,t}.E(ii));
    
    qdd = [M11,M12;M21,M22] \([check.shoulder{1,s,t}.torque_c(ii); check.elbow{1,s,t}.torque_c(ii)] - [C1;C2]);
    
%     qdd(1) = I \(check.shoulder{1,s,t}.torque_c(ii) - C1);
%     qdd(1) = I \(check.elbow{1,s,t}.torque_c(ii) - C2);
    
    check.theta{1,s,t}.Sdd(ii) = qdd(1);
    check.theta{1,s,t}.Edd(ii) = qdd(2);
end
% vars{1,s,t}.masses
%% Set Init Conditions
check.theta{1,s,t}.S(1,:) = theta{1,s,t}.S(1,:);
check.theta{1,s,t}.E(1,:) = theta{1,s,t}.E(1,:);
check.theta{1,s,t}.Sd(1,1:length(theta{1,s,t}.Sd(1,:))) = theta{1,s,t}.Sd(1,:);
check.theta{1,s,t}.Ed(1,1:length(theta{1,s,t}.Ed(1,:))) = theta{1,s,t}.Ed(1,:);

%% Calulate Kinematics
check.theta{1,s,t}.Sd = (cumsum(check.theta{1,s,t}.Sdd)*vars{1,s,t}.time_inc+check.theta{1,s,t}.Sd(1))';
check.theta{1,s,t}.Ed = (cumsum(check.theta{1,s,t}.Edd)*vars{1,s,t}.time_inc+check.theta{1,s,t}.Ed(1))';

check.theta{1,s,t}.S = (cumsum(check.theta{1,s,t}.Sd)*vars{1,s,t}.time_inc+check.theta{1,s,t}.S(1));
check.theta{1,s,t}.E = (cumsum(check.theta{1,s,t}.Ed)*vars{1,s,t}.time_inc+check.theta{1,s,t}.E(1));

check.theta{1,s,t}.Sdd = check.theta{1,s,t}.Sdd';
check.theta{1,s,t}.Edd = check.theta{1,s,t}.Edd';

check.theta{1,s,t}.S = reshape(check.theta{1,s,t}.S,[length(check.theta{1,s,t}.S),1]);
check.theta{1,s,t}.E = reshape(check.theta{1,s,t}.E,[length(check.theta{1,s,t}.E),1]);

check.theta{1,s,t}.Sd = reshape(check.theta{1,s,t}.Sd,[length(check.theta{1,s,t}.Sd),1]);
check.theta{1,s,t}.Ed = reshape(check.theta{1,s,t}.Ed,[length(check.theta{1,s,t}.Ed),1]);

check.theta{1,s,t}.Sdd = reshape(check.theta{1,s,t}.Sdd,[length(check.theta{1,s,t}.Sdd),1]);
check.theta{1,s,t}.Edd = reshape(check.theta{1,s,t}.Edd,[length(check.theta{1,s,t}.Edd),1]);

check.x{1,s,t} = upperarm{1,s,t}.length*cos(check.theta{1,s,t}.S)+...
    forearm{1,s,t}.length*cos(check.theta{1,s,t}.S+check.theta{1,s,t}.E);
check.y{1,s,t} = sin(check.theta{1,s,t}.S)*upperarm{1,s,t}.length+...
    sin(check.theta{1,s,t}.E+check.theta{1,s,t}.S)*forearm{1,s,t}.length;

check.vx{1,s,t} = diff(check.x{1,s,t})/.0025;
check.vy{1,s,t} = diff(check.y{1,s,t})/.0025;
v1(loop,31:length(check.vx{1,s,t}.^2+check.vy{1,s,t}.^2)+30)...
    = sqrt(check.vx{1,s,t}.^2+check.vy{1,s,t}.^2);
switch c
    case 1
        c1(loop) = c;
        s1(loop) = s;
        t1(loop) = t;
    case 2
        c1(loop) = c;
        s1(loop) = s;
        t1(loop) = t;
    case 3
        c1(loop) = c;
        s1(loop) = s-1;
        t1(loop) = t; 
    case 4
        c1(loop) = c;
        s1(loop) = s-1;
        t1(loop) = t;
end
% check.v(1,s,t) = sqrt(check.vx{1,s,t}.^2+check.vy{1,s,t}.^2)

%% Calc Error

s_error(1,s,t,d) = sum(abs((check.theta{1,s,t}.S - theta{1,s,t}.S)));
e_error(1,s,t,d) = sum(abs((check.theta{1,s,t}.E - theta{1,s,t}.E)));
figure(2);hold on
plot(v1(loop,:))
% drawnow;
[c,s,t]
loop = loop + 1

        end
    end
end


for s=2:6
    for t = 1:4
        for c=3:4
for ii = 1:length(act{1,s-1,t}.an)
    norm_force = vars{1,s-1,t}.norm_force;
    for k = 1:length(muscle_nums)
        norm_length(k) = muscles{1,s-1,t}.(muscle_nums{k}).length(ii)/muscles{1,s-1,t}.(muscle_nums{k}).l0;
        vel(k) = muscles{1,s-1,t}.(muscle_nums{k}).v(ii);
        switch c
            case 1
                a = act{1,s-1,t}.(muscle_nums{k})(ii);
            case 2
                a = act{1,s-1,t}.(muscle_nums{k})(ii);
            case 3
                a = act{1,s-1,t}.(muscle_nums{k})(ii);
            case 4
                a = act{1,s-1,t}.(muscle_nums{k})(ii);                
        end
        
        check.stress(k,ii) = Fl_Fv_for(norm_length(k),vel(k),a)*norm_force;
        check.force(k,ii) = check.stress(k,ii)*muscles{1,s-1,t}.(muscle_nums{k}).pcsa;
    end
    
    A1 = [muscles{1,s-1,t}.an.m_arm_e(ii),...
        muscles{1,s-1,t}.bs.m_arm_e(ii),...
        muscles{1,s-1,t}.br.m_arm_e(ii),...
        0,...
        0,...
        0,...
        muscles{1,s-1,t}.bb.m_arm_e(ii),...
        muscles{1,s-1,t}.tb.m_arm_e(ii)];
    A2 = [0,...
        0,...
        0,...
        muscles{1,s-1,t}.da.m_arm_s(ii),...
        muscles{1,s-1,t}.dp.m_arm_s(ii),...
        muscles{1,s-1,t}.pc.m_arm_s(ii),...
        muscles{1,s-1,t}.bb.m_arm_s(ii),...
        muscles{1,s-1,t}.tb.m_arm_s(ii)];
    A=[A1;A2];
    
    x = [check.force(1,ii),check.force(2,ii),check.force(3,ii),...
        check.force(4,ii),check.force(5,ii),check.force(6,ii),...
        check.force(7,ii),check.force(8,ii)]';

    check.elbow{1,s-1,t}.torque_c(ii) = A1*x;
    check.shoulder{1,s-1,t}.torque_c(ii) = A2*x;
    
    prm.m1 = upperarm{1,s-1,t}.mass;
    prm.r1 = upperarm{1,s-1,t}.centl;
    prm.l1 = upperarm{1,s-1,t}.length;
    prm.i1 = upperarm{1,s-1,t}.Ic;
    
    prm.m2 = forearm{1,s-1,t}.mass;
    prm.r2 = forearm{1,s-1,t}.l_com;
    prm.r22 = forearm{1,s-1,t}.centl;
    prm.l2 = forearm{1,s-1,t}.length;
    prm.i2 = forearm{1,s-1,t}.Ic;
    
    prm.m = vars{1,s-1,t}.masses;
    
    M11 = prm.m1*prm.r1^2 + prm.i1 +...
        (prm.m+prm.m2)*(prm.l1^2+prm.r22^2+...
        2*prm.l1*prm.r22*cos(theta{1,s-1,t}.E(ii))) +...
        prm.i2;

    M12 = (prm.m2+prm.m)*(prm.r22^2+...
        prm.l1*prm.r22*cos(theta{1,s-1,t}.E(ii))) +...
        prm.i2;

    M21 = M12;

    M22 = prm.m2*prm.r2^2+prm.m*prm.l2^2+prm.i2;
    
    C1 = -forearm{1,s-1,t}.mass*forearm{1,s-1,t}.centl*upperarm{1,s-1,t}.length*(theta{1,s-1,t}.Ed(ii).^2).*sin(theta{1,s-1,t}.E(ii))-...
        2*forearm{1,s-1,t}.mass*forearm{1,s-1,t}.centl*upperarm{1,s-1,t}.length*theta{1,s-1,t}.Sd(ii).*theta{1,s-1,t}.Ed(ii).*sin(theta{1,s-1,t}.E(ii));

    C2 = forearm{1,s-1,t}.mass*forearm{1,s-1,t}.centl*upperarm{1,s-1,t}.length*(theta{1,s-1,t}.Sd(ii).^2).*sin(theta{1,s-1,t}.E(ii));
    
    qdd = [M11,M12;M21,M22] \([check.shoulder{1,s-1,t}.torque_c(ii); check.elbow{1,s-1,t}.torque_c(ii)] - [C1;C2]);
    
%     qdd(1) = I \(check.shoulder{1,s-1,t}.torque_c(ii) - C1);
%     qdd(1) = I \(check.elbow{1,s-1,t}.torque_c(ii) - C2);
    
    check.theta{1,s-1,t}.Sdd(ii) = qdd(1);
    check.theta{1,s-1,t}.Edd(ii) = qdd(2);
end
% vars{1,s-1,t}.masses
%% Set Init Conditions
check.theta{1,s-1,t}.S(1,:) = theta{1,s-1,t}.S(1,:);
check.theta{1,s-1,t}.E(1,:) = theta{1,s-1,t}.E(1,:);
check.theta{1,s-1,t}.Sd(1,1:length(theta{1,s-1,t}.Sd(1,:))) = theta{1,s-1,t}.Sd(1,:);
check.theta{1,s-1,t}.Ed(1,1:length(theta{1,s-1,t}.Ed(1,:))) = theta{1,s-1,t}.Ed(1,:);

%% Calulate Kinematics
check.theta{1,s-1,t}.Sd = (cumsum(check.theta{1,s-1,t}.Sdd)*vars{1,s-1,t}.time_inc+check.theta{1,s-1,t}.Sd(1))';
check.theta{1,s-1,t}.Ed = (cumsum(check.theta{1,s-1,t}.Edd)*vars{1,s-1,t}.time_inc+check.theta{1,s-1,t}.Ed(1))';

check.theta{1,s-1,t}.S = (cumsum(check.theta{1,s-1,t}.Sd)*vars{1,s-1,t}.time_inc+check.theta{1,s-1,t}.S(1));
check.theta{1,s-1,t}.E = (cumsum(check.theta{1,s-1,t}.Ed)*vars{1,s-1,t}.time_inc+check.theta{1,s-1,t}.E(1));

check.theta{1,s-1,t}.Sdd = check.theta{1,s-1,t}.Sdd';
check.theta{1,s-1,t}.Edd = check.theta{1,s-1,t}.Edd';

check.theta{1,s-1,t}.S = reshape(check.theta{1,s-1,t}.S,[length(check.theta{1,s-1,t}.S),1]);
check.theta{1,s-1,t}.E = reshape(check.theta{1,s-1,t}.E,[length(check.theta{1,s-1,t}.E),1]);

check.theta{1,s-1,t}.Sd = reshape(check.theta{1,s-1,t}.Sd,[length(check.theta{1,s-1,t}.Sd),1]);
check.theta{1,s-1,t}.Ed = reshape(check.theta{1,s-1,t}.Ed,[length(check.theta{1,s-1,t}.Ed),1]);

check.theta{1,s-1,t}.Sdd = reshape(check.theta{1,s-1,t}.Sdd,[length(check.theta{1,s-1,t}.Sdd),1]);
check.theta{1,s-1,t}.Edd = reshape(check.theta{1,s-1,t}.Edd,[length(check.theta{1,s-1,t}.Edd),1]);

check.x{1,s-1,t} = upperarm{1,s-1,t}.length*cos(check.theta{1,s-1,t}.S)+...
    forearm{1,s-1,t}.length*cos(check.theta{1,s-1,t}.S+check.theta{1,s-1,t}.E);
check.y{1,s-1,t} = sin(check.theta{1,s-1,t}.S)*upperarm{1,s-1,t}.length+...
    sin(check.theta{1,s-1,t}.E+check.theta{1,s-1,t}.S)*forearm{1,s-1,t}.length;

check.vx{1,s-1,t} = diff(check.x{1,s-1,t})/.0025;
check.vy{1,s-1,t} = diff(check.y{1,s-1,t})/.0025;
v1(loop,31:length(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2)+30)...
    = sqrt(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2);
switch c
    case 1
        c1(loop) = c;
        s1(loop) = s;
        t1(loop) = t;
    case 2
        c1(loop) = c;
        s1(loop) = s;
        t1(loop) = t;
    case 3
        c1(loop) = c;
        s1(loop) = s-1;
        t1(loop) = t; 
    case 4
        c1(loop) = c;
        s1(loop) = s-1;
        t1(loop) = t;
end
% check.v(1,s-1,t) = sqrt(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2)

%% Calc Error

s_error(1,s-1,t,d) = sum(abs((check.theta{1,s-1,t}.S - theta{1,s-1,t}.S)));
e_error(1,s-1,t,d) = sum(abs((check.theta{1,s-1,t}.E - theta{1,s-1,t}.E)));
figure(2);hold on
plot(v1(loop,:))
% drawnow;
[c,s,t]
loop = loop + 1

        end
    end
end

cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Data');
fileID=fopen('react_sim.csv','w');
A = {'c,','s,','t,','method,','react_indx'};
fprintf(fileID,'%s',A{:});
fprintf(fileID,'\n');

for k = 1:size(v1,1)
    v = v1(k,:);
%     v_diff = diff23f5(v,1/200,12);
    v_diff = diff(v);
    ind1 = find(v>.03,1,'first');
    % std(Data.TanV(j-10:j,i))<6e-4 || (std(V_diff(j-10:j,i))<1e-4 && std(Data.TanV(j-10:j,i))<2e-3)
    while ind1 > 11
        if std(v(ind1-10:ind1))<2e-3 
            break
        end
        ind1 = ind1 - 1;
    end
    react_indx1(k) = ind1;
    
    ind2 = find(v>.03,1,'first');
    while ind2 > 11
        if std(v(ind2-10:ind2))<6e-4 || (std(v_diff(ind2-10:ind2))<1e-4 && std(v(ind2-10:ind2))<2e-3)
            break
        end
        ind2 = ind2 - 2;
    end
    react_indx2(k) = ind2;
    
    A = {c1(k),s1(k),t1(k),3,react_indx1(k)};
    dlmwrite('react_sim.csv',A,'-append');
    A = {c1(k),s1(k),t1(k),4,react_indx2(k)};
    dlmwrite('react_sim.csv',A,'-append');
end

fclose(fileID);


