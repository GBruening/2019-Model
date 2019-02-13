cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files');
clear a_hold a b fval
% if ~exist('est')
%     load('metspeed_73117_highnormforce_cap.mat');
% end
% if ~exist('usum1')
%     Plotting;
% end
% Plotting;

if ~exist('titlestr')
    titlestr = filename;
    graphstr = filename;
    titlestr2 = filename;
    graphstr2 = filename;
end

% plot3 = 1;
% if plot3
%     figure(3);clf(3);
% end
d=2;
conf_int = .975;
ColorSet2=parula(18);

func2 = @(a,x,m,t) (a(1)*(m^a(2))*(.1)^(a(3)))/(t^a(4));
init_guess = [1,1,1,1];
a_min = [0,0,0,0];
a_max = [1E5,1E5,1E5,1E5];

[num,txt,raw] = xlsread('Met_Data.xlsx');
[a,b] = size(num);
for ii=1:a
%     txt{ii+2,1}
    if ~isempty(char(txt{ii+2,1})) && length(char(txt{ii+2,1}))>3
        if strcmp(char(txt{ii+2,1}(2:4)),'mp0')
            c=1;
            for s = 1:6
                met_data(c,s) = num(ii-1,s);
            end
        end
        if strcmp(char(txt{ii+2,1}(2:4)),'mp5')
            c=2;
            for s = 1:6
                met_data(c,s) = num(ii-1,s);
            end
        end
        if strcmp(char(txt{ii+2,1}(2:4)),'mp1')
            c=3;
            for s = 1:6
                met_data(c,s) = num(ii-1,s);
            end
        end
        if strcmp(char(txt{ii+2,1}(2:4)),'mp2')
            c=4;
            for s = 1:6
                met_data(c,s) = num(ii-1,s);
            end
        end
    end
end

% opts = optimset('Display','off');
% x=ones(size(met_data));
% func4 = @(a,x) func3;
% [a,~,resid,~,output,~,J] = ...
%     lsqcurvefit(func4,init_guess,ones(size(met_data)),met_data,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(bsum2),max(bsum2),100),a,resid,'Jacobian',J);
% 
% bounds = nlparci(a,resid,'Jacobian',J);


func2 = @(a,x) func_fit2(a,x,vars,met_data);
opts = optimset('Display','off');
[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,[1,1,1,1,1],ones(size(met_data)),met_data,[0,0,0,0,0],[100,100,100,100,100],opts);

y_pred = func4(a,1);

figure(5);clf(5); hold on
for c = 1:4
    for s = 1:6
        ColorPref(c) = plot(vars{c,s,t}.speeds(s),met_data(c,s),...
            'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        plot(vars{c,s,t}.speeds(s),y_pred(c,s),...
            'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
    end
    plot(vars{c,s,t}.speeds,met_data(c,:),'k-');
    plot(vars{c,s,t}.speeds,y_pred(c,:),'--','Color',ColorSet(c,:));
end
y = met_data;
y_fit = y_pred;
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;
    
text(.9, .8*max(max(y_pred)), sprintf('%0.3f+%0.3f*(m^{%0.3f})*(d^{%0.3f})/(T^{%0.3f}) \nR^2 = %0.2f'...
    ,a(1),a(2),a(3),a(4),a(5),R_sq));
xlabel('Time');ylabel('Umberger Model');
    legend([ColorPref(:)],mass_legend);

% func3 = @(a,x) func4(a,x);
% [y_pred,delta] = nlpredci(func3,bsum1,a,resid,'Jacobian',J);
% bounds = nlparci(a,resid,'Jacobian',J);
% fprintf('IntFit  a=%.2f, lb = %.2f, ub = %.2f\n',a,bounds(1),bounds(2));

function [y_pred] = func_fit2(a,x,vars,y)
    
    for m = 1:4
        switch m
            case 1
            t = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
            mass = 2.470;
            case 2
            t = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
            mass = 4.8550;
            case 3
            t = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
            mass = 7.1600;
            case 4
            t = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
            mass = 11.7250;
        end    
        y_pred(m,:) = a(1) + (a(2).*(mass.^a(3)).*(.1).^(a(4)))./(t.^a(5));
    end
    
    error = sum(sum((y-y_pred).^2));
end

function [y_pred] = func4(a,x)
    
    for m = 1:4
        switch m
            case 1
                t = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
                mass = 2.470;
            case 2
                t = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
                mass = 4.8550;
            case 3
                t = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
                mass = 7.1600;
            case 4
                t = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
                mass = 11.7250;
        end    
        y_pred(m,:) = a(1) + (a(2).*(mass.^a(3)).*(.1).^(a(4)))./(t.^a(5));
    end
%     y_pred= y_pred(:)
%     error = sum(sum((y-y_pred).^2));
end









