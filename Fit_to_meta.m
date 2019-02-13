% cd('E:\Documents\OneDrive\Muscle modeling\Min_jerk files');
cd(main_folder)
clear a_hold a b fval
% if ~exist('est')
%     cd('D:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files\Data');
%     load('force2cost_21-Nov-2017.mat');
%     cd('D:\Users\Gary\OneDrive\Muscle modeling\Min_jerk files');
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

plot3 = 1;
if plot3
    figure(3);clf(3);
end
d=2;
conf_int = .975;
ColorSet2=parula(18);

if ~exist('fit_method')
    fit_method = 'free';
    p = 1;
    expo = 1;
    k1 = 7;
end

switch fit_method
    case 'linear'
        expo_min = 1;
        expo_max = 1;
        func2 = @(a,x) a(1)*x.^1+a(2);
        init_guess = [1;1];
        a_min = [-1E5,-1E5];
        a_max = [1E5,1E5];
        func3 = @(a,x) func_fit_linear;
    case 'squared'
        expo_min = 2;
        expo_max = 2;
        func2 = @(a,x) a(1)*x.^2+a(2);
        init_guess = [1;1];
        a_min = [-1E5,-1E5];
        a_max = [1E5,1E5];
        func3 = @(a,x) func_fit_squared;
    case 'free'
        expo_min = 0;
        expo_max = 10;
        func2 = @(a,x) a(1)*x.^a(2)+a(3);     
        init_guess = [1;1;1];  
        a_min = [-1E5,-1E5,-1E5];
        a_max = [1E5,1E5,1E5];
        func3 = @(a,x) func_fit_free;
end

% [num,txt,raw] = xlsread('Met_Data.xlsx');
% [a,b] = size(num);
% for ii=1:a
% %     txt{ii+2,1}
%     if ~isempty(char(txt{ii+2,1})) && length(char(txt{ii+2,1}))>3
%         if strcmp(char(txt{ii+2,1}(2:4)),'mp0')
%             c=1;
%             for s = 1:6
%                 met_data(c,s) = num(ii-1,s);
%                 met_data_err(c,s) = num(ii,s);
%             end
%         end
%         if strcmp(char(txt{ii+2,1}(2:4)),'mp5')
%             c=2;
%             for s = 1:6
%                 met_data(c,s) = num(ii-1,s);
%                 met_data_err(c,s) = num(ii,s);
%             end
%         end
%         if strcmp(char(txt{ii+2,1}(2:4)),'mp1')
%             c=3;
%             for s = 1:6
%                 met_data(c,s) = num(ii-1,s);
%                 met_data_err(c,s) = num(ii,s);
%             end
%         end
%         if strcmp(char(txt{ii+2,1}(2:4)),'mp2')
%             c=4;
%             for s = 1:6
%                 met_data(c,s) = num(ii-1,s);
%                 met_data_err(c,s) = num(ii,s);
%             end
%         end
%     end
% end

% met_data = csvread('met_power_data.csv',1,0);
mp_data = readtable('met_power_data.csv');


% metpowernet
% test = met_data.movedur(met_data.subject==1 & met_data.eff_mass==2.47)

eff_masses = [2.47,4.73,6.99,11.5];
met_data = zeros([4,8,7]);
for c=1:4
    for subj = 1:8
        tempdurs = mp_data.movedur(mp_data.subject==subj & mp_data.effmass==eff_masses(c))
        tempmet = mp_data.metpowernet(mp_data.subject==subj & mp_data.effmass==eff_masses(c));
        for k = 1:length(tempdurs)
            counter = counter+1;
            [~,s]=min(abs(tempdurs(k)-vars{c,subj,3,t}.speeds));
            if isempty(vars{c,subj,s,t})
                continue
            end            
            met_data(c,subj,s) = tempmet(k);
        end
    end
end

figure(2);clf(2);
%% Neural Drive
figure(2);
h=subplot(2,2,1);
ii=1;
clear x y
clear a_hold a
for c = 1:4
    for s = 1:6
        
        if ~isnan(met_data(c,s))
        
            usum2(ii) = usum1(c,s);
            y(ii) = met_data(c,s);
            ii=ii+1;
        end
    end
end
opts = optimset('Display','off');
[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,init_guess,usum2,y,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(usum2),max(usum2),100),a,resid,'Jacobian',J);

% min_fit = y_pred-delta;
% max_fit = y_pred+delta;

ci  = nlparci(a,resid,'Jacobian',J);

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

x_pred = linspace(min(usum2),max(usum2),100);
y_fit = a(1).*usum2.^a(2)+a(3);
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;

RSQ_matrix(p,expo,k1,1) = R_sq;

% func3 = @(a,x) func_fit3(a,x);
nonlinfit = fitnlm(usum2,y,func2,init_guess);
R_sq = nonlinfit.Rsquared.Ordinary;
for k = 1:length(init_guess)
    a(k) = nonlinfit.Coefficients{k,1};
    se(k) = nonlinfit.Coefficients{k,2};
end

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end
RSQ_matrix(p,expo,k1,1) = R_sq;

fit = a(1).*usum1.^a(2) + a(3);
y_pred = a(1).*x_pred.^a(2) + a(3);
% min_fit = ci(1).*x_pred.^ci(2)+ci(3);
% max_fit = ci(4).*x_pred.^ci(5)+ci(6);

if exist('plotting','var') && strcmp(plotting,'sensitivity')
    figure(5);
    h=subplot(2,5,count); hold on;
else
    figure(2);
    subplot(2,2,1);hold on;
end
hold on
% shadedErrorBar(x_pred,y_pred,[max_fit;min_fit],...
%     {'Color',ColorSet2(15,:),'Color',ColorSet2(15,:)},.3);
plot(x_pred,y_pred,'k-');
plot(usum2,y,'x','Color',ColorSet2(6,:));
eq_str = sprintf('%.2f*x^{%.2f}+%.2f = y \nr^{2} = %.2f'...
    ,a(1),a(2),a(3),R_sq);
xlabel(sprintf('Neural Drive Rate^%g',expo));
ylabel('Metabolic Cost???');
title(titlestr);
text(h.XLim(1)+(h.XLim(2)-h.XLim(1))/10,h.YLim(2)*.8,eq_str);

if plot3
    figure(3);
    subplot(2,2,1);hold on;
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
            plot(vars{c,subj,s,t}.speeds(s),fit(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        errorbar(vars{c,subj,s,t}.speeds,met_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         plot(vars{c,subj,s,t}.speeds,met_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,fit(c,:),'--','Color',ColorSet(c,:));
    end
         pbaspect([1 1 1]);     xlim([.2 1.4]);

    legend([ColorPref(:)],mass_legend);
    xlabel('Movement Duration');ylabel(sprintf('Cost (drive^%g predict)',expo)); 
    title(titlestr2);
end

%% Active State
figure(2);
h=subplot(2,2,2);
ii=1;
if length(a_max)== 3
    a_max(3) = 40;
else    
    a_max(2) = 40;
end
clear x y
clear a_hold a
for c = 1:4
    for s = 1:6
        
        if ~isnan(met_data(c,s))
        
            asum2(ii) = asum1(c,s);
            y(ii) = met_data(c,s);
            ii=ii+1;
        end
    end
end

opts = optimset('Display','off');
[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,init_guess,asum2,y,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(asum2),max(asum2),100),a,resid,'Jacobian',J);

% min_fit = y_pred-delta;
% max_fit = y_pred+delta;

ci  = nlparci(a,resid,'Jacobian',J);

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

x_pred = linspace(min(asum2),max(asum2),100);
y_fit = a(1).*asum2.^a(2)+a(3);
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;
RSQ_matrix(p,expo,k1,2) = R_sq;

if expo ==2;
    1;
end
func3 = @(a,x) func_fit3(a,x);
nonlinfit = fitnlm(asum2,y,func2,init_guess);
R_sq = nonlinfit.Rsquared.Ordinary;
for k = 1:length(init_guess)
    a(k) = nonlinfit.Coefficients{k,1};
    se(k) = nonlinfit.Coefficients{k,2};
end

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end
RSQ_matrix(p,expo,k1,2) = R_sq;

fit = a(1).*asum1.^a(2) + a(3);
y_pred = a(1).*x_pred.^a(2) + a(3);
% min_fit = ci(1).*x_pred.^ci(2)+ci(3);
% max_fit = ci(4).*x_pred.^ci(5)+ci(6);

hold on
if exist('plotting','var') && strcmp(plotting,'sensitivity')
    figure(6);
    h=subplot(2,5,count); hold on;
else
    figure(2);
    subplot(2,2,2);hold on;
end
plot(x_pred,y_pred,'k-');
plot(asum2,y,'x','Color',ColorSet2(6,:));
eq_str = sprintf('%.2f*x^{%.2f}+%.2f = y \nr^{2} = %.2f'...
    ,a(1),a(2),a(3),R_sq);
xlabel(sprintf('Active State Rate^%g',expo));
ylabel('Metabolic Cost???');
title(titlestr);
text(h.XLim(1)+(h.XLim(2)-h.XLim(1))/10,h.YLim(2)*.8,eq_str);

if plot3
    
    
    figure(3);
    subplot(2,2,2);hold on;
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
            plot(vars{c,subj,s,t}.speeds(s),fit(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        errorbar(vars{c,subj,s,t}.speeds,met_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         plot(vars{c,subj,s,t}.speeds,met_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,fit(c,:),'--','Color',ColorSet(c,:));
    end
         pbaspect([1 1 1]);     xlim([.2 1.4]);

    legend([ColorPref(:)],mass_legend);
    xlabel('Movement Duration');ylabel(sprintf('Cost (Act State^%g Predict)',expo)); 
    title(titlestr2);
    
end

%% Muscle Force
figure(2);
h=subplot(2,2,3);
ii=1;
if length(a_max)== 3
    a_max(3) = 1E5;
else    
    a_max(2) = 1E5;
end
clear x y
clear a_hold a
for c = 1:4
    for s = 1:6
        
        if ~isnan(met_data(c,s))
        
            fsum2(ii) = fsum1(c,s);
            y(ii) = met_data(c,s);
            ii=ii+1;
        end
    end
end

1;
opts = optimset('Display','off');
% [beta,R,J] = nlinfit(fsum2,y,@func2,init_guess);

% fsum3 = fsum2/1E4;

[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,init_guess,fsum2,y,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(fsum3),max(fsum3),100),a,resid,'Jacobian',J);

% min_fit = y_pred-delta;
% max_fit = y_pred+delta;

ci  = nlparci(a,resid,'Jacobian',J);


switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

x_pred = linspace(min(fsum2),max(fsum2),100);
y_fit = a(1).*fsum2.^a(2)+a(3);
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;
RSQ_matrix(p,expo,k1,3) = R_sq;

% func3 = @(a,x) func_fit3(a,x);
nonlinfit = fitnlm(fsum2,y,func2,init_guess);
R_sq = nonlinfit.Rsquared.Ordinary;
for k = 1:length(init_guess)
    a(k) = nonlinfit.Coefficients{k,1};
    se(k) = nonlinfit.Coefficients{k,2};
end

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end
RSQ_matrix(p,expo,k1,3) = R_sq;

fit = a(1).*(fsum1).^a(2) + a(3);
y_pred = a(1).*x_pred.^a(2) + a(3);
% min_fit = ci(1).*x_pred.^ci(2)+ci(3);
% max_fit = ci(4).*x_pred.^ci(5)+ci(6);

hold on
if exist('plotting','var') && strcmp(plotting,'sensitivity');
    figure(7);
    h=subplot(2,5,count); hold on;
else
    figure(2);
    subplot(2,2,3);hold on;
end
plot(x_pred,y_pred,'k-');
plot(fsum2,y,'x','Color',ColorSet2(6,:));

eq_str = sprintf('%.2fx^{%.2f}+%.2f = y \nr^{2} = %.2f'...
    ,a(1),a(2),a(3),R_sq);
xlabel(sprintf('Muscle Force Rate^%g',expo));
ylabel('Metabolic Cost???');
title(titlestr);
text(h.XLim(1)+(h.XLim(2)-h.XLim(1))/10,h.YLim(2)*.8,eq_str);

if plot3
    fit = a(1).*(fsum1).^a(2) + a(3);
    
    figure(3);
    subplot(2,2,3);hold on;
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
            plot(vars{c,subj,s,t}.speeds(s),fit(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        errorbar(vars{c,subj,s,t}.speeds,met_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         plot(vars{c,subj,s,t}.speeds,met_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,fit(c,:),'--','Color',ColorSet(c,:));
    end
         pbaspect([1 1 1]);     xlim([.2 1.4]);

    legend([ColorPref(:)],mass_legend);
    xlabel('Movement Duration');ylabel(sprintf('Cost (Force Rate^%g Predict)',expo)); 
    title(titlestr2);
    
end

%% Joint Torque
figure(2);
h=subplot(2,2,4);
ii=1;
clear x y
clear a_hold a
for c = 1:4
    for s = 1:6
        
        if ~isnan(met_data(c,s))
        
            tsum2(ii) = tsum1(c,s);
            y(ii) = met_data(c,s);
            ii=ii+1;
        end
    end
end


opts = optimset('Display','off');
[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,init_guess,tsum2,y,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(tsum2),max(tsum2),100),a,resid,'Jacobian',J);
% 
% min_fit = y_pred-delta;
% max_fit = y_pred+delta;

ci  = nlparci(a,resid,'Jacobian',J);


switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

x_pred = linspace(min(tsum2),max(tsum2),100);
y_fit = a(1).*tsum2.^a(2)+a(3);
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;
RSQ_matrix(p,expo,k1,4) = R_sq;

% func3 = @(a,x) func_fit3(a,x);
nonlinfit = fitnlm(tsum2,y,func2,init_guess);
R_sq = nonlinfit.Rsquared.Ordinary;
for k = 1:length(init_guess)
    a(k) = nonlinfit.Coefficients{k,1};
    se(k) = nonlinfit.Coefficients{k,2};
end

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end
RSQ_matrix(p,expo,k1,4) = R_sq;

fit = a(1).*tsum1.^a(2) + a(3);
y_pred = a(1).*x_pred.^a(2) + a(3);
% min_fit = ci(1).*x_pred.^ci(2)+ci(3);
% max_fit = ci(4).*x_pred.^ci(5)+ci(6);

hold on
if exist('plotting','var') && strcmp(plotting,'sensitivity');
    figure(8);
    h=subplot(2,5,count); hold on;
else
    figure(2);
    subplot(2,2,4);hold on;
end
plot(x_pred,y_pred,'k-');
plot(tsum2,y,'x','Color',ColorSet2(6,:));
eq_str = sprintf('%.2f*x^{%.2f}+%.2f = y \nr^{2} = %.2f'...
    ,a(1),a(2),a(3),R_sq);
xlabel(sprintf('Joint Torque Rate^%g',expo));
ylabel('Metabolic Cost???');
title(titlestr);
text(h.XLim(1)+(h.XLim(2)-h.XLim(1))/10,h.YLim(2)*.8,eq_str);

fit = a(1).*tsum1.^a(2) + a(3);
sqrt(sum(sum((met_data-fit).^2)/(4*6-2)))/(sqrt(sum((usum2-mean(usum2)).^2)));

if plot3
    
    fit = a(1).*tsum1.^a(2) + a(3);
    
    figure(3);
    subplot(2,2,4);hold on;
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
            plot(vars{c,subj,s,t}.speeds(s),fit(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        errorbar(vars{c,subj,s,t}.speeds,met_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         plot(vars{c,subj,s,t}.speeds,met_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,fit(c,:),'--','Color',ColorSet(c,:));
    end
         pbaspect([1 1 1]);     xlim([.2 1.4]);

    legend([ColorPref(:)],mass_legend);
    xlabel('Movement Duration');ylabel(sprintf('Cost (Joint Torque^%g predict)',expo)); 
    title(titlestr2);
    
end


%% Umberger model
figure(45);clf(45);
h=subplot(1,1,1);
ii=1;
clear x y
clear a_hold a
for c = 1:4
    for s = 1:6
        
        if ~isnan(met_data(c,s))
        
            bsum2(ii) = bsum1(c,s);
            y(ii) = met_data(c,s);
            ii=ii+1;
        end
    end
end
opts = optimset('Display','off');
[a,~,resid,~,output,~,J] = ...
    lsqcurvefit(func2,init_guess,bsum2,y,a_min,a_max,opts);
% [y_pred,delta] = nlpredci(func2,linspace(min(bsum2),max(bsum2),100),a,resid,'Jacobian',J);
% 
% min_fit = y_pred-delta;
% max_fit = y_pred+delta;

% ci  = nlparci(a,resid,'Jacobian',J);

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

x_pred = linspace(min(bsum2),max(bsum2),100);
y_fit = a(1).*(bsum2).^a(2)+a(3);
SSR = sum((y_fit-mean(y)).^2);
SSTO = sum((y-mean(y)).^2);
R_sq = SSR/SSTO;
RSQ_matrix(p,expo,k1,5) = R_sq;


% func3 = @(a,x) func_fit3(a,x);
nonlinfit = fitnlm(bsum2,y,func2,init_guess);
R_sq = nonlinfit.Rsquared.Ordinary;
for k = 1:length(init_guess)
    a(k) = nonlinfit.Coefficients{k,1};
    se(k) = nonlinfit.Coefficients{k,2};
end

switch fit_method
    case 'linear'        
        a(3) = a(2);
        a(2) = 1;
    case 'squared'
        a(3) = a(2);
        a(2) = 2;
end

fit = a(1).*(bsum1).^a(2) + a(3);
y_pred = a(1).*x_pred.^a(2) + a(3);
% min_fit = ci(1).*x_pred.^ci(2)+ci(3);
% max_fit = ci(4).*x_pred.^ci(5)+ci(6);

if exist('plotting','var') && strcmp(plotting,'sensitivity');
    figure(9);
    h=subplot(2,5,count); hold on;
else
    figure(45);
    h=subplot(2,1,1);hold on;
end
hold on
plot(x_pred,y_pred,'k-');
plot(bsum2,y,'x','Color',ColorSet2(6,:));
eq_str = sprintf('%.2f*x^{%.2f}+%.2f = y \nr^{2} = %.2f'...
    ,a(1),a(2),a(3),R_sq);
xlabel(sprintf('Umberger^%g',expo));
ylabel('Metabolic Cost???');
title(titlestr);
text(h.XLim(1)+(h.XLim(2)-h.XLim(1))/10,h.YLim(2)*.8,eq_str);
         pbaspect([1 1 1]);

if plot3
    figure(45);
    subplot(2,1,2);
    hold on;
%     subplot(2,2,1);hold on;
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:),'linewidth',1);
            plot(vars{c,subj,s,t}.speeds(s),fit(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        errorbar(vars{c,subj,s,t}.speeds,met_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         plot(vars{c,subj,s,t}.speeds,met_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,fit(c,:),'--','Color',ColorSet(c,:));
    end
         pbaspect([1 1 1]);     xlim([.2 1.4]);

    legend([ColorPref(:)],mass_legend);
    xlabel('Movement Duration');ylabel(sprintf('Cost')); 
    title(titlestr2);
    
end

% figure(56);
% clf(56);
% hold on
% for c = 1:4
%     for s = 1:6
%     ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),met_data(c,s),...
%             'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
%     end
%     plot(vars{c,subj,s,t}.speeds,met_data(c,:),'Color',ColorSet(c,:));
% end
% legend([ColorPref(:)],mass_legend);
% xlabel('Movement Duration');ylabel(sprintf('Cost (Metabolics)')); 
% title(titlestr2);
    
% figure(57);
% clf(57);
% hold on
% for c = 1:4
%     for s = 1:6
%     ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),bsum1(c,s),...
%             'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
%     end
%     plot(vars{c,subj,s,t}.speeds,bsum1(c,:),'--','Color',ColorSet(c,:));
% end
% legend([ColorPref(:)],mass_legend);
% xlabel('Movement Duration');ylabel(sprintf('Cost (Umberger)')); 
% title(titlestr2);
%     


%% New fitting

% func2 = @(a,x,m,t) (a(1)*(m^a(2))*(.1)^(a(3)))/(t^a(4));
init_guess = [1,1,1,1];
a_min = [0,0,0,0];
a_max = [1E5,1E5,1E5,1E5];

% if p == 1 && expo == 1

ii=1;
proxies = {'tsum1','fsum1','asum1','usum1','bsum1','met_data'};
proxy_labels = {'Torque','Force','Active State','Neural Drive','Umberger','Metabolic Power'};
figure(4);clf(4); hold on
for v = 1:length(proxies)
    clear x y
    clear a_hold a
    switch proxies{v}
        case 'tsum1'
            fit_data = tsum1;
        case 'fsum1'
            fit_data = fsum1;
        case 'asum1'
            fit_data = asum1;
        case 'usum1'
            fit_data = usum1;
        case 'bsum1'
            fit_data = bsum1;
        case 'met_data'
            fit_data = met_data;
            1;
    end
    func2 = @(a,x) func_fit2(a,x,vars,met_data);
    opts = optimset('Display','off');
    [a,~,resid,~,output,~,J] = ...
        lsqcurvefit(func2,[1,1,1,1,1],ones(size(met_data)),fit_data,[0,0,0,0,0],[100,100,100,100,100],opts);
    ci  = nlparci(a,resid,'Jacobian',J);
    y_pred = func4(a,vars);

    1;
    func3 = @(a,x) func_fit3(a,x);
    for c = 1:4
        tum(c,:) = vars{c,s,1}.speeds;
    end
    for s = 1:6
        muss(:,s) = [2.470; 4.8550; 7.1600; 11.7250];
    end
    x = [tum(:),muss(:)];
    try
        nonlinfit = fitnlm(x,fit_data(:),func3,[20,2,1,.7,3.2]);
    catch
        nonlinfit = fitnlm(x,fit_data(:),func3,10*rand(5,1));
    end
    R_sq = nonlinfit.Rsquared.Ordinary;
    for k = 1:5
        a(k) = nonlinfit.Coefficients{k,1};
        se(k) = nonlinfit.Coefficients{k,2};
    end
    y_pred = func4(a,vars);
%     clear time
%     [y_pred,time] = func_manytime(a,vars);
    
%     y_fit = h.
%     y = fit_data;
%     y_fit = y_pred;
%     SSR = sum((y_fit-mean(y)).^2);
%     SSTO = sum((y-mean(y)).^2);
%     R_sq = SSR/SSTO;
    Comp_matrix{k1,v}.Rsq = R_sq;
    Comp_matrix{k1,v}.a = a;
    Comp_matrix{k1,v}.conf_int = ci;
    Comp_matrix{k1,v}.se = se;

    subplot(2,3,v);
    hold on
    for c = 1:4
        for s = 1:6
            ColorPref(c) = plot(vars{c,subj,s,t}.speeds(s),fit_data(c,s),...
                'o','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
            plot(vars{c,subj,s,t}.speeds(s),y_pred(c,s),...
                'x','Color',ColorSet(c,:),'MarkerFaceColor',ColorSet(c,:));
        end
        plot(vars{c,subj,s,t}.speeds,fit_data(c,:),'k-');
        plot(vars{c,subj,s,t}.speeds,y_pred(c,:),'Color',ColorSet(c,:));
%         plot(time(c,:),time(c,:).*y_pred(c,:),'Color',ColorSet(c,:));
%         if v == 6
%             errorbar(vars{c,subj,s,t}.speeds,vars{c,subj,s,t}.speeds.*fit_data(c,:),met_data_err(c,:)./sqrt(8),'Color',ColorSet(c,:));
%         end
    end
    text(.9, .8*max(max(y_pred)), ...
        sprintf('%s = %0.3f+%0.3f*(m^{%0.3f}))/(T^{%0.3f}) \n R^2 = %0.2f',...
        proxy_labels{v},a(1),a(2),a(3),a(4),R_sq));
    xlabel('Time');ylabel(proxy_labels{v});
        legend([ColorPref(:)],mass_legend);
%          pbaspect([1 1 1]); 
         xlim([.2 1.4]);
end
cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Graphs\Compare Fits');
figure(4);
beautifyfig;
drawnow;
% % savefig(minparams{k1});
% end
% 
x_data = [tsum2',fsum2',asum2',fsum2',bsum2'];
% y_data = [y,met_data_err(:)];


1;
%% Saving Figs
% drawnow;
% cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files\Graphs\Cost Fits');
% if ~exist('plotting','var')
%     figure(2);
%     beautifyfig;
%     savefig(graphstr);
%     if exist(graphstr,'file')==2
%         fprintf('Graph saved as %s \n',graphstr');
%     end
%     cd([main_folder filesep 'Graphs' filesep 'Money Plots']);
%     figure(3);
%     beautifyfig;
%     savefig(graphstr2);
%     figure(45);
%     beautifyfig;
%     cd([main_folder filesep 'Graphs' filesep 'Umberfits']);
%     savefig(umbstr);
% else
%     cd([main_folder filesep 'Graphs' filesep 'Sensitivity']);
%     figure(5);
%     savefig('Drive');
%     figure(6);
%     savefig('State');
%     figure(7);
%     savefig('force');
%     figure(8);
%     savefig('torque');
%     figure(9);
%     savefig('Umberger');
% %     savefig(graphstr);
%     fprintf('Graph saved as: sensitivty %s \n',graphstr');
% end
% 
cd(main_folder);
% 1;
% % 
% % %% Write files to text doc
cd('D:\Users\Gary\Google Drive\Muscle modeling\Min_jerk files');
% fileID = fopen('Summed Data.txt');
% 
% for c = 1:4
%     for s= 1:6
%         fprintf(fileID,'%g %g %g %g %g %f %f %f %f %f %f %f %f %f %f %f \n',...
%             p,expo,k1,vars{c,s,1}.masses,s,vars{c,s,1}.speeds(s),...
%             tsum1(c,s),RSQ_matrix(p,expo,k1,1),fsum1(c,s),RSQ_matrix(p,expo,k1,2),...
%             asum1(c,s),RSQ_matrix(p,expo,k1,3),usum1(c,s),RSQ_matrix(p,expo,k1,4),...
%             bsum1(c,s),RSQ_matrix(p,expo,k1,5));
%     end
% end

%% Mix
% figure(3);
% % subplot(2,2,4);
% ii=1;
% clear x y
% clear a_hold a
% for c = 1:4
%     for s = 1:6
%         
%         if ~isnan(met_data(c,s))
%         
%             mix1(ii) = usum1(c,s);
%             mix2(ii) = tsum1(c,s);
%             mix2(ii) = 0;
%             y(ii) = met_data(c,s);
%             ii=ii+1;
%         end
%     end
% end
% 
% func = @(a) func_fit_mix(a,mix1,mix2,y);
% options = optimset('TolFun',1E-8,'TolX',1E-8,'Display','off');
% for run = 1:20
%     [a_hold(run,:),fval(run),exitflag,output] = ...
%         fmincon(func,[100*rand(1)*ones(1,5)],...
%         [],[],[],[],[0,0,0,0,0],[1E5,1E5,1E5,1E5,1E5],[], options);
% end
% [a1,b]=min(fval);
% a=a_hold(b,:);
% 
% x1_pred = linspace(min(mix1),max(mix1),100);
% x2_pred = linspace(min(mix2),max(mix2),100);
% y_pred = a(1).*x1_pred.^a(2)+a(3).*x2_pred.^a(4)+a(5);
% y_fit = a(1).*mix1.^a(2)+a(3).*mix2.^a(4)+a(5);
% R_sq = 1-min(fval)

%%
function [y_pred] = func_fit2(a,x,vars,y)
    
    for m = 1:4
        t = vars{m,1,1}.speeds;
        switch m
            case 1
%             t = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
            mass = 2.470;
            case 2
%             t = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
            mass = 4.8550;
            case 3
%             t = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
            mass = 7.1600;
            case 4
%             t = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
            mass = 11.7250;
        end    
        y_pred(m,:) = a(1) + (a(2).*(mass.^a(3)))./(t.^a(4));
    end
    
    error = sum(sum((y-y_pred).^2));
end

function [y_pred] = func_fit3(a,x)
    m = x(:,2);
    t = x(:,1);
    y_pred = a(1) + (a(2).*(.1^a(3)).*(m.^a(4)))./(t.^a(5));
end

function [y_pred] = func_fit_linear(a,x)
    y_pred = a(1).*x.^1+a(2);
end

function [y_pred] = func_fit_squared(a,x)
    y_pred = a(1).*x.^2+a(2);
end

function [y_pred] = func_fit_free(a,x)
    y_pred = a(1).*x.^a(2)+a(3);
end

function [y_pred] = func4(a,vars)
    
    for m = 1:4
        t = vars{m,1,1}.speeds;
        switch m
            case 1
%                 t = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
                mass = 2.470;
            case 2
%                 t = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
                mass = 4.8550;
            case 3
%                 t = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
                mass = 7.1600;
            case 4
%                 t = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
                mass = 11.7250;
        end    
        y_pred(m,:) = a(1) + (a(2).*(.1^a(3)).*(mass.^a(4)))./(t.^a(5));
    end
%     y_pred= y_pred(:)
%     error = sum(sum((y-y_pred).^2));
end

function [y_pred,time] = func_manytime(a,vars)
    
    for m = 1:4
        switch m
            case 1
%                 t = [0.4407, 0.492, 0.5908, 0.7733, 0.965, 1.1567];
                time(m,:) = linspace(vars{m,1,1}.speeds(1)-.05,1.5,100);
                mass = 2.470;
            case 2
%                 t = [0.4703,	0.5089,	0.5956,	0.7807,	0.9694,	1.1657];
                time(m,:) = linspace(vars{m,1,1}.speeds(1)-.05,1.5,100);
                mass = 4.8550;
            case 3
%                 t = [0.5114,	0.6006,	0.7847,	0.9701,	1.1591,	1.3438];
                time(m,:) = linspace(vars{m,1,1}.speeds(1)-.05,1.5,100);
                mass = 7.1600;
            case 4
%                 t = [0.5301,	0.6049,	0.7849,	0.9756,	1.1634,	1.3371,];
                time(m,:) = linspace(vars{m,1,1}.speeds(1)-.05,1.5,100);
                mass = 11.7250;
        end    
        y_pred(m,:) = a(1) + (a(2).*(mass.^a(3)))./(time(m,:).^a(4));
    end
%     y_pred= y_pred(:)
%     error = sum(sum((y-y_pred).^2));
end

function [error] = func_fit(a,x,y)
    
    func = @(a,x) a(1).*x.^a(2) + a(3);
    
    y_pred = func(a,x);
    
    error = sum((y-y_pred).^2);
    
end


function [error] = func_fit_mix(a,x1,x2,y)
    
    func = @(a,x1,x2) a(1).*x1.^a(2) + a(3).*x2.^a(4) + a(5);
    
    y_pred = func(a,x1,x2);
    
    error = sum((y-y_pred).^2);
    
    x1_pred = linspace(min(x1),max(x1),100);
    x2_pred = linspace(min(x2),max(x2),100);
    y_pred = a(1).*x1_pred.^a(2)+a(3).*x2_pred.^a(4)+a(5);
    y_fit = a(1).*x1.^a(2)+a(3).*x2.^a(4)+a(5);
    SSR = sum((y_fit-mean(y)).^2);
    SSTO = sum((y-mean(y)).^2);
    SSE = sum((y-y_fit).^2);
%     R_sq1 = 1-SSE/SSTO;
    
    R_sq = 1 - sum((y - y_fit).^2)/sum((y - mean(y)).^2);
    
%     R_sq = 1-SSR/SSTO;
    
    error = 1-R_sq;
    
end