
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
v1(loop,1:length(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2))...
    = sqrt(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2);
switch c
    case 1
        c1(loop) = c;
        s2(loop) = s;
        t1(loop) = t;
    case 2
        c1(loop) = c;
        s2(loop) = s;
        t1(loop) = t;
    case 3
        c1(loop) = c;
        s2(loop) = s-1;
        t1(loop) = t; 
    case 4
        c1(loop) = c;
        s2(loop) = s-1;
        t1(loop) = t;
end
% check.v(1,s-1,t) = sqrt(check.vx{1,s-1,t}.^2+check.vy{1,s-1,t}.^2)

%% Calc Error

s_error(1,s-1,t,d) = sum(abs((check.theta{1,s-1,t}.S - theta{1,s-1,t}.S)));
e_error(1,s-1,t,d) = sum(abs((check.theta{1,s-1,t}.E - theta{1,s-1,t}.E)));
figure(2);hold on
plot(v(loop,:))
drawnow;
[c,s,t]
loop = loop + 1

        end
    end
end