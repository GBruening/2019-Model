function [ shoulder , elbow , theta , eff_mass] = ...
    Calc_kine_test( Data , forearm , upperarm, time_inc , add_mass)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

order = 5;
window = 21;

shoulder.torque=zeros(max(size(Data.x,1)),size(Data.x,2));
elbow.torque=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.E=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.S=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.Sd=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.Ed=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.Sdd=zeros(max(size(Data.x,1)),size(Data.x,2));
theta.Edd=zeros(max(size(Data.x,1)),size(Data.x,2));

trial = 1;

% for trial=1:1%size(Data.x,2)
    
%% Get joint angles
index = ~isnan(Data.x(:));

x=Data.x(~isnan(Data.x(index)));
y=Data.y(~isnan(Data.y(index)));

a = x.^2 + y.^2 - upperarm.length.^2 - forearm.length.^2;
b = 2 * upperarm.length * forearm.length;

theta.E(index) = atan2(sqrt(1-(a.^2./b)),a./b);

k1 = upperarm.length + forearm.length*cos(theta.E(index));
k2 = forearm.length*sin(theta.E(index));
theta.S(index) = atan2(y,x)-atan2(k2,k1);

s1 = diff23f5(theta.S,.0025,3);
e1 = diff23f5(theta.E,.0025,3);
theta.Sd = s1(:,2);
theta.Ed = e1(:,2);
theta.Sd(2) = (theta.Sd(1)+theta.Sd(3))/2;
theta.Ed(2) = (theta.Ed(1)+theta.Ed(3))/2;

esp = 16;
theta.Sd = [interp1([1 2 esp-2:esp'], [0 0 theta.Sd(esp-2:esp)'], [1:esp], 'spline')';...
    theta.Sd(esp+1:end)];
theta.Ed = [interp1([1 2 esp-2:esp'], [0 0 theta.Ed(esp-2:esp)'], [1:esp], 'spline')';...
    theta.Ed(esp+1:end)];
theta.Sdd = [s1(:,3)];
theta.Edd = [e1(:,3)];
theta.Sdd(2) = (theta.Sdd(1)+theta.Sdd(3))/2;
theta.Edd(2) = (theta.Edd(1)+theta.Edd(3))/2;
theta.Sdd = [interp1([1 2 esp-2:esp'], [0 0 theta.Sdd(esp-2:esp)'], [1:esp], 'spline')';...
    theta.Sdd(esp+1:end)];
theta.Edd = [interp1([1 2 esp-2:esp'], [0 0 theta.Edd(esp-2:esp)'], [1:esp], 'spline')';...
    theta.Edd(esp+1:end)];

a1 = (upperarm.mass * (upperarm.centl)^2 + upperarm.Ic);
a2 = (forearm.mass * (forearm.centl)^2 + forearm.Ic);
a3 = (forearm.mass * upperarm.length * abs(forearm.centl));
a4 = (forearm.mass * upperarm.length^2);

% Alaa's q1 is the shoulder angle
% Joint torques

prm.m1 = upperarm.mass;
prm.r1 = upperarm.centl;
prm.l1 = upperarm.length;
prm.i1 = upperarm.Ic;

prm.m2 = forearm.mass;
prm.r2 = forearm.l_com;
prm.r22 = forearm.centl;
prm.l2 = forearm.length;
prm.i2 = forearm.Ic;

prm.m = add_mass;

M11 = prm.m1*prm.r1^2 + prm.i1 +...
    (prm.m+prm.m2)*(prm.l1^2+prm.r22^2+...
    2*prm.l1*prm.r22*cos(theta.E(:))) +...
    prm.i2;

M12 = (prm.m2+prm.m)*(prm.r22^2+...
    prm.l1*prm.r22*cos(theta.E(:))) +...
    prm.i2;

M21 = M12;

M22 = prm.m2*prm.r2^2+prm.m*prm.l2^2+prm.i2;

C1 = -forearm.mass*forearm.centl*upperarm.length*(theta.Ed(:).^2).*sin(theta.E(:))-...
    2*forearm.mass*forearm.centl*upperarm.length*theta.Sd(:).*theta.Ed(:).*sin(theta.E(:));

C2 = forearm.mass*forearm.centl*upperarm.length*(theta.Sd(:).^2).*sin(theta.E(:));

T1 = M11.*theta.Sdd(:)+M12.*theta.Edd(:)+C1;
T2 = M21.*theta.Sdd(:)+M22.*theta.Edd(:)+C2;

shoulder.torque_a(:) = T1;
elbow.torque_a(:) = T2;

shoulder.torque(:) = shoulder.torque_a(:);
elbow.torque(:) = elbow.torque_a(:);

for ii = 1:length(M11)
    J11 = -upperarm.length.*sin(theta.S(ii))-forearm.length*sin(theta.S(ii)+theta.E(ii));
    J12 = -forearm.length*sin(theta.S(ii)+theta.E(ii));
    J21 = upperarm.length.*cos(theta.S(ii))+forearm.length*cos(theta.S(ii)+theta.E(ii));
    J22 = forearm.length*cos(theta.S(ii)+theta.E(ii));

    J = [J11 J12 ; J21 J22];
    arm.J = J;
    arm.I=[M11(ii) M12(ii); M21(ii) M22];
    arm.M=transpose(inv(arm.J))*arm.I*inv(arm.J);
    move_angle = atan2(Data.targetposition(1)-Data.startposition(1),...
        Data.targetposition(2)-Data.startposition(2));
    eff_mass(ii) = norm(arm.M*[cos(move_angle) sin(move_angle)]')+1;
end
    
  1;
    %% uncomment to line 248 for checks
%      x1=upperarm.length*cos(theta.S(:)) + ...
%         forearm.length*cos(theta.S(:)+theta.E(:));
%     y1=upperarm.length*sin(theta.S(:)) + ...
%         forearm.length*sin(theta.S(:)+theta.E(:));
% 
%     x1=x1(y1>0);
%     y1=y1(y1>0);
% 
%     figure(1);clf(1);
%     hold on
%     plot(Data.x(:),Data.y(:),'o','Color','r')
%     plot(x1,y1,'*','Color','g');
% %     viscircles(ro,vars.distances(d),'color','k');
%     plot(Data.targetposition(1),Data.targetposition(2),'x','Linewidth',5);
%     1;
%     %% Check torque's lining up with calc'd angles
%     figure(1);figure(2);
%     clf(1);clf(2);
%     for k=1:length(elbow.torque(:))
%         clear I C qd T
%         I(1,1) = a1+a4+2*a3*cos(theta.E(k))+a2;
%         I(1,2) = a2+a3*cos(theta.E(k));
%         I(2,1) = a2+a3*cos(theta.E(k));
%         I(2,2) = a2;
%         
%         C(1,1) = -a3*theta.Ed(k)*sin(theta.E(k));
%         C(1,2) = -a3*(theta.Ed(k)+theta.Sd(k))*sin(theta.E(k));
%         C(2,1) = a3*sin(theta.E(k))*theta.Sd(k);
%         C(2,2) = 0;
%         
%         qd(1,1) = theta.Sd(k);
%         qd(2,1) = theta.Ed(k);
%         
%         T(1,1) = shoulder.torque(k);
%         T(2,1) = elbow.torque(k);
%         
%         qdd = I \ (T-C*qd);
%         
% %         figure(1);hold on
% %         plot(k,qdd(1),'x','Color','b');
% %         
% %         figure(2);hold on
% %         plot(k,qdd(2),'x','Color','b');
%         
%         theta.Sdd_c(k) = qdd(1)';
%         theta.Edd_c(k) = qdd(2)';
%         
%     end
%     
%     figure(1); hold on
%     plot(theta.Sdd_c,'x');plot(theta.Sdd(:,1)); 
%     legend({'Test','inv kin'});
%     figure(2);hold on
%     plot(theta.Edd_c,'x');plot(theta.Edd(:,1)); legend({'Test','inv kin'});
%     
%     check.theta.S(1,:) = theta.S(1,:);
%     check.theta.E(1,:) = theta.E(1,:);
%     check.theta.Sd(1,1:length(theta.Sd(1,:))) = theta.Sd(1,:);
%     check.theta.Ed(1,1:length(theta.Ed(1,:))) = theta.Ed(1,:);
% 
%     %% Calulate Kinematics
%     check.theta.Sdd = theta.Sdd_c;
%     check.theta.Edd = theta.Edd_c;
%     check.theta.Sd = (cumsum(check.theta.Sdd)*.0025+check.theta.Sd(1))';
%     check.theta.Ed = (cumsum(check.theta.Edd)*.0025+check.theta.Ed(1))';
% 
%     check.theta.S = (cumsum(check.theta.Sd)*.0025+check.theta.S(1));
%     check.theta.E = (cumsum(check.theta.Ed)*.0025+check.theta.E(1));
% 
%     check.theta.Sdd = check.theta.Sdd';
%     check.theta.Edd = check.theta.Edd';
% 
%     check.theta.S = reshape(check.theta.S,[length(check.theta.S),1]);
%     check.theta.E = reshape(check.theta.E,[length(check.theta.E),1]);
% 
%     check.theta.Sd = reshape(check.theta.Sd,[length(check.theta.Sd),1]);
%     check.theta.Ed = reshape(check.theta.Ed,[length(check.theta.Ed),1]);
% 
%     check.theta.Sdd = reshape(check.theta.Sdd,[length(check.theta.Sdd),1]);
%     check.theta.Edd = reshape(check.theta.Edd,[length(check.theta.Edd),1]);
%     figure(2);clf(2);
% hold on
% plot(check.theta.S)
% plot(theta.S)
%     1;
%     if abs(theta.S(end)-(sum(cumsum(theta.Sdd_c))*time_inc^2+theta.S(1)))>1E-7 ||...
%             abs(theta.E(end)-(sum(cumsum(theta.Edd_c))*time_inc^2+theta.E(1)))>1E-7
%     
%         1;
%     end
    
    %% Alaa Method for Angle
%     for ii=1:length(index)
%         
%         dist(ii) = sqrt(Data.x(ii)^2+Data.y(ii)^2);
%         
%         theta.S_a(ii) = (acos((Data.x(ii)^2-Data.y(ii)^2-...
%             upperarm.length^2-forearm.length^2)/(2*upperarm.length*forearm.length)));
%         
%         theta.b_a(ii) = asin(forearm.length*sin(theta.S(ii))/dist(ii));
%         theta.c_a(ii) = asin(Data.x(ii)/dist(ii));
%         
%         if Data.x(ii)<0
%             theta.E_a(ii) = pi - theta.c_a(ii) - theta.b_a(ii);
%         else
%             theta.E_a(ii) = theta.c_a(ii) - theta.b_a(ii);
%         end
%         
%     end
%     
%     % Jacobian For joint kinematics


%     Hxd = sgolayfilt(Data.vx(:),order,window);
%     Hyd = sgolayfilt(Data.vy(:),order,window);
    
%     Hx = [ones(1,floor(window/2))'*Hxd(1);Hxd];
%     Hy = [ones(1,floor(window/2))'*Hyd(1);Hyd];

%     Hx = Data.vx(:);
%     Hy = Data.vy(:);
%     clear Hxd Hxy
%     
%     for ii = 1:sum(index)
%         
%         if ii == 260
%             1;
%         end
%         
%         J(1,1) = -upperarm.length.*sin(theta.S(ii))-forearm.length*sin(theta.S(ii)+theta.E(ii));
%         J(1,2) = -forearm.length*sin(theta.S(ii)+theta.E(ii));
%         J(2,1) = upperarm.length.*cos(theta.S(ii))+forearm.length*cos(theta.S(ii)+theta.E(ii));
%         J(2,2) = forearm.length*cos(theta.S(ii)+theta.E(ii));
%         
%         
%         E = J\[Hx(ii);Hy(ii)];
%         
%         theta.Sd_J(ii) = E(1) * .005;
%         theta.Ed_J(ii) = E(2) * .005;
%         
%         Jd(1,1)= -upperarm.length*theta.Sd(ii)*cos(theta.S(ii))-...
%             forearm.length*(theta.Sd(ii)+theta.Ed(ii))*cos(theta.S(ii)+theta.E(ii));
%         Jd(1,2)= -forearm.length*(theta.Sd(ii)+theta.Ed(ii))*cos(theta.S(ii)+theta.E(ii));
%         Jd(2,1)= -upperarm.length*theta.Sd(ii)*sin(theta.S(ii))-...
%             forearm.length*(theta.Sd(ii)+theta.Ed(ii))*sin(theta.S(ii)+theta.E(ii));
%         Jd(2,2)= -forearm.length*(theta.Sd(ii)+theta.Ed(ii))*sin(theta.S(ii)+theta.E(ii));
%         
%         
%         
%         clear E
% %         E = inv(Jd) * ([a1;
%     end
    
    
    
    %%
%     a = [a1+a4,a3.*cos(theta.E(index)-theta.S(index));a3.*cos(theta.E(index)-theta.S(index)),a2]';
%     b = [theta.Sdd(index)';theta.Edd(index)'];
%     b = repmat(b,length(x),1);
%     c = [zeros(length(x),1),-a3.*theta.Ed(index).*sin(theta.E(index)-theta.S(index));a3.*theta.Sd(index).*sin(theta.E(index)-theta.S(index)),a2]';
%     d = [theta.Sd(index)';theta.Ed(index)'];
%     d = repmat(d,length(x),1);
%     e = a*b + c*d;
% 
%     shoulder.torque(1:length(e(1,:))) = e(1,:);
%     elbow.torque(1:length(e(1,:))) = e(2,:);

% end


end

