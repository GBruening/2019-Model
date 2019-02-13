function[time, x, y, vx, vy, ax, ay]=Gen_mvt_gb(Resamp,subj,c,s,t, ro, rf, time_inc)
1;
% if c == 3 || c == 4
%     s1 = s+1;
% else
%     s1 = s;
% end

s1=s;

x = reshape(((Resamp.P(c,s1,:)-Resamp.P(c,s1,1))/max(Resamp.P(c,s1,:)))*(rf(1)-ro(1)),100,1)+ro(1);
y = reshape(((Resamp.P(c,s1,:)-Resamp.P(c,s1,1))/max(Resamp.P(c,s1,:)))*(rf(2)-ro(2)),100,1)+ro(2);
% x1 = diff23f5(x,.0025,50);
% x1(1,2) = 0;
% x1(1,3) = 0;
% y1 = diff23f5(y,.0025,50);
% y1(1,2) = 0;
% y1(1,3) = 0;
% plot(sqrt((x1(:,1)-ro(1)).^2+(y1(:,1)-ro(2)).^2))

buff = 0.005;
x1 = spline([0,.0025,reshape(buff+Resamp.T(c,s1,4:end),1,length(Resamp.T(c,s1,4:end)))],...
    [ro(1);ro(1);x(4:end)],[0:time_inc:Resamp.T(c,s1,end)])';
y1 = spline([0,.0025,reshape(buff+Resamp.T(c,s1,4:end),1,length(Resamp.T(c,s1,4:end)))],...
    [ro(2);ro(2);y(4:end)],[0:time_inc:Resamp.T(c,s1,end)])';
x1(1) = ro(1);
y1(1) = ro(2);

% x1 = sgolayfilt(x1,1,9);
% y1 = sgolayfilt(x1,1,9);

x1 = diff23f5(x1,.0025,10);
x1(1,2) = 0;
x1(1,3) = 0;
y1 = diff23f5(y1,.0025,10);
y1(1,2) = 0;
y1(1,3) = 0;

x = reshape(x1(:,1),length(x1(:,1)),1);
y = reshape(y1(:,1),length(y1(:,1)),1);

x = [x(1);x(1);x(1);x];
y = [y(1);y(1);y(1);y];

vx = x1(:,2)*.0025;
vy = y1(:,2)*.0025;
ax = x1(:,3)*.0025*.0025;
ay = y1(:,3)*.0025*.0025;
% vx=0;
% vy=0;
% ax=0;
% ay=0;


% plot(time,sqrt((x1-ro(1)).^2+(y1-ro(2)).^2))

time = time_inc*[0:length(x)];%[0:time_inc:Resamp.T(c,s1,end)+.0075];
time = reshape(time,length(time),1);
% c
% s
% 
% plot(time,sqrt((x1(:,1)-ro(1)).^2+(y1(:,1)-ro(2)).^2))
% 1;