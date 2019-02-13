theta.E = linspace(0,150,150);
theta.S = ones(size(theta.E))*45;
% theta.S = linspace(0,90,90);

muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};

Colors = [...
         0  ,  0.4470  ,  0.7410;...
    0.8500  ,  0.3250  ,  0.0980;...
    0.9290  ,  0.6940  ,  0.1250;...
    0.4940  ,  0.1840  ,  0.5560;...
    0.4660  ,  0.6740  ,  0.1880;...
    0.3010  ,  0.7450  ,  0.9330;...
    0.6350  ,  0.0780  ,  0.1840;...
         0  ,  0.4470  ,  .3000];

%from langenderfer 2004
muscles.an.pcsa = 1.89E-4;
muscles.bs.pcsa = 8.71E-4;  
muscles.br.pcsa = 2.15E-4;  
muscles.da.pcsa = 6.46E-4;
muscles.dp.pcsa = 5.69E-4;
muscles.pc.pcsa = 6.68E-4;
muscles.bb.pcsa = (1.57+1.75)*1E-4;
muscles.tb.pcsa = (4.13+3.60+3.21)*1E-4;

%FT from https://ac.els-cdn.com/0022510X73900233/1-s2.0-0022510X73900233-main.pdf?_tid=f2587a6a-d3c7-11e7-b5e6-00000aacb360&acdnat=1511824317_c65b8060cce54a277f6130de38be9eac
% or http://www.sciencedirect.com/science/article/pii/S0021929004005147?via%3Dihub
muscles.an.ft = .4;
muscles.bs.ft = .4;
muscles.br.ft = .602;
muscles.da.ft = .467/2+.390/2;
muscles.dp.ft = .467/2+.390/2;
muscles.pc.ft = .57;
muscles.bb.ft = .577/2+.495/2;
muscles.tb.ft = .59/2+.467/2;

muscles.an.coef_e = [-2.7306E-9,10.448E-7,-14.329E-5,8.4297E-3,-2.2841E-1,-5.3450]; 
muscles.bs.coef_e = [-2.0530E-5,2.3425E-3,2.3080E-1,5.5492];
muscles.br.coef_e = [-5.6171E-5,10.084E-3,1.6681E-1,19.49];
muscles.da.coef_s = 33.02;
muscles.dp.coef_s = -78.74;
muscles.pc.coef_s = 50.80;
muscles.bb.coef_s = 29.21;
muscles.bb.coef_e = [-2.9883E-5,1.8047E-3,4.5322E-1,14.660];
muscles.tb.coef_s = -25.40;
muscles.tb.coef_e = [-3.5171E-9,13.277E-7,-19.092E-5,12.886E-3,-3.0284E-1,-23.287];


muscles.an.m_arm_e = polyval(muscles.an.coef_e,theta.E)/1000; %Divide by 10
muscles.bs.m_arm_e = polyval(muscles.bs.coef_e,theta.E)/1000; %covnert mm to
muscles.br.m_arm_e = polyval(muscles.br.coef_e,theta.E)/1000; %cm
muscles.da.m_arm_s = polyval(muscles.da.coef_s,theta.S)/1000;
muscles.dp.m_arm_s = polyval(muscles.dp.coef_s,theta.S)/1000;
muscles.pc.m_arm_s = polyval(muscles.pc.coef_s,theta.S)/1000;
muscles.bb.m_arm_s = polyval(muscles.bb.coef_s,theta.S)/1000;
muscles.bb.m_arm_e = polyval(muscles.bb.coef_e,theta.E)/1000;
muscles.tb.m_arm_s = polyval(muscles.tb.coef_s,theta.S)/1000;
muscles.tb.m_arm_e = polyval(muscles.tb.coef_e,theta.E)/1000;

figure(1);
clf(1)
subplot(2,1,1);hold on
plot(muscles.an.m_arm_e,'Color',Colors(1,:),'linewidth',3);
plot(muscles.bs.m_arm_e,'Color',Colors(2,:),'linewidth',3);
plot(muscles.br.m_arm_e,'Color',Colors(3,:),'linewidth',3);
plot(muscles.bb.m_arm_e,'Color',Colors(7,:),'linewidth',3);
plot(muscles.tb.m_arm_e,'Color',Colors(8,:),'linewidth',3);
legend({'an','bs','br','bb','tb'});
xlabel('Elbow Angle (deg)');
set(gca,'ytick',[-0.04:.02:.08],'yticklabel',[-4:2:8])
ylabel('Moment arm (cm)')
subplot(2,1,2);hold on
plot(muscles.da.m_arm_s,'Color',Colors(4,:),'linewidth',3);
plot(muscles.dp.m_arm_s,'Color',Colors(5,:),'linewidth',3);
plot(muscles.pc.m_arm_s,'Color',Colors(6,:),'linewidth',3);
plot(muscles.bb.m_arm_s,'Color',Colors(7,:),'linewidth',3);
plot(muscles.tb.m_arm_s,'Color',Colors(8,:),'linewidth',3);
legend({'da','dp','pc','bb','tb'});
xlabel('Shoulder Angle (deg)');
set(gca,'ytick',[-0.1:.02:.08],'yticklabel',[-10:2:8])
ylim([-.1 .06])
ylabel('Moment arm (cm)')

muscles.an.coef_l = [4.7658E-11,-1.8235E-8,25.008E-7,-14.713E-5,3.9865E-3,9.3288E-2,0];
muscles.bs.coef_l = [3.5832E-7,-4.0884E-5,-4.0282E-3,-9.6852E-2,0];
muscles.br.coef_l = [11.374E-7,-17.60E-5,-2.9114E-3,-34.017E-2,0];
muscles.da.coef_l = [-5.7631E-1,0];
muscles.dp.coef_l = [13.743E-1,0];
muscles.pc.coef_l = [-8.8663E-1,0];
muscles.bb.coef_l = [5.2156E-7,-3.1498E-5,-7.9101E-3,-25.587E-2,0];
muscles.tb.coef_l = [6.1385E-11,-2.3174E-8,33.321E-7,-22.491E-5,5.2856E-3,40.644E-2,0];

%Divide muscle lengths by 1000 to convert mm to m
muscles.an.length = (polyval(muscles.an.coef_l,theta.E) + 53.57)/1000;
muscles.bs.length = (polyval(muscles.bs.coef_l,theta.E) + 137.48)/1000;
muscles.br.length = (polyval(muscles.br.coef_l,theta.E) + 276.13)/1000;
muscles.da.length = (polyval(muscles.da.coef_l,theta.S) + 172.84)/1000;
muscles.dp.length = (polyval(muscles.dp.coef_l,theta.S) + 157.64)/1000;
muscles.pc.length = (polyval(muscles.pc.coef_l,theta.S) + 155.19)/1000;
muscles.bb.length = (polyval(muscles.bb.coef_l,theta.E) + -5.0981E-1*theta.S+378.06)/1000;
muscles.tb.length = (polyval(muscles.tb.coef_l,theta.E) + 4.4331E-1*theta.S +260.05)/1000;

figure(1);
clf(1)
hold on
for k = 1:length(muscle_nums)
%  h(k) = 
 plot(muscles.(muscle_nums{k}).length,'linewidth',3,'color',Colors(k,:));
%  h(k).Color
end
xlabel('Elbow Angle (deg)');
ylabel('Muscle Length (m)');
title('$\theta_S = 45^\circ$','Interpreter','latex');
legend(muscle_nums);