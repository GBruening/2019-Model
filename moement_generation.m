vars.norm_force = 220e4
%from langenderfer 2004
% muscles1.an.pcsa = 2.89E-4;
% muscles1.bs.pcsa = 8.71E-4;  
% muscles1.br.pcsa = 2.15E-4;  
% muscles1.da.pcsa = 6.46E-4;
% muscles1.dp.pcsa = 5.69E-4;
% muscles1.pc.pcsa = 6.68E-4;
% muscles1.bb.pcsa = (1.57+1.75)*1E-4;
% muscles1.tb.pcsa = (4.13+3.60+3.21)*1E-4;

muscles1.an.pcsa = 1.3E-4;
muscles1.bs.pcsa = 7.71E-4;
muscles1.br.pcsa = 1.15E-4;  
muscles1.da.pcsa = 6.46E-4;
muscles1.dp.pcsa = 3.00E-4;
% muscles1.pc.pcsa = 6.07E-4;
muscles1.pc.pcsa = 6.07E-4;
muscles1.bb.pcsa = (1.57+1.75)*1E-4;
muscles1.tb.pcsa = (4.13+3.60+3.21-1)*1E-4;

% muscles1.an.pcsa = 1.3E-4;
% muscles1.bs.pcsa = 7.71E-4;  
% muscles1.br.pcsa = 1.15E-4;  
% muscles1.da.pcsa = 5.46E-4;
% muscles1.dp.pcsa = 4.69E-4;
% muscles1.pc.pcsa = 3.07E-4;
% 
% muscles1.bb.pcsa = (1.57+1.75)*1E-4;
% muscles1.tb.pcsa = (4.13+3.60+3.21)*1E-4;

muscles1.an.coef_e = [-2.7306E-9,10.448E-7,-14.329E-5,8.4297E-3,-2.2841E-1,-5.3450]; 
muscles1.bs.coef_e = [-2.0530E-5,2.3425E-3,2.3080E-1,5.5492];
muscles1.br.coef_e = [-5.6171E-5,10.084E-3,1.6681E-1,19.49];
muscles1.da.coef_s = 33.02;
muscles1.dp.coef_s = -78.74;
muscles1.pc.coef_s = 50.80;
muscles1.bb.coef_s = 29.21;
muscles1.bb.coef_e = [-2.9883E-5,1.8047E-3,4.5322E-1,14.660];
muscles1.tb.coef_s = -25.40;
muscles1.tb.coef_e = [-3.5171E-9,13.277E-7,-19.092E-5,12.886E-3,-3.0284E-1,-23.287];

muscles1.an.coef_e = [-2.7306E-9,10.448E-7,-14.329E-5,8.4297E-3,-2.2841E-1,-5.3450]; 
muscles1.bs.coef_e = [-2.0530E-5,2.3425E-3,2.3080E-1,5.5492];
muscles1.br.coef_e = [-5.6171E-5,10.084E-3,1.6681E-1,19.49];
muscles1.da.coef_s = 33.02;
muscles1.dp.coef_s = -78.74;
muscles1.pc.coef_s = 50.80;
muscles1.bb.coef_s = 29.21;
muscles1.bb.coef_e = [-2.9883E-5,1.8047E-3,4.5322E-1,14.660];
muscles1.tb.coef_s = -25.40;
muscles1.tb.coef_e = [-3.5171E-9,13.277E-7,-19.092E-5,12.886E-3,-3.0284E-1,-23.287];


muscles1.an.m_arm_e = 2*polyval(muscles1.an.coef_e,90)/1000; %Divide by 10
muscles1.bs.m_arm_e = polyval(muscles1.bs.coef_e,90)/1000; %covnert mm to
muscles1.br.m_arm_e = polyval(muscles1.br.coef_e,90)/1000; %cm
muscles1.da.m_arm_s = polyval(muscles1.da.coef_s,45)/1000;
muscles1.dp.m_arm_s = polyval(muscles1.dp.coef_s,45)/1000;
muscles1.pc.m_arm_s = polyval(muscles1.pc.coef_s,45)/1000;
muscles1.bb.m_arm_s = polyval(muscles1.bb.coef_s,45)/1000;
muscles1.bb.m_arm_e = polyval(muscles1.bb.coef_e,90)/1000;
muscles1.tb.m_arm_s = polyval(muscles1.tb.coef_s,45)/1000;
muscles1.tb.m_arm_e = 2*polyval(muscles1.tb.coef_e,90)/1000;

n_t = Fl_Fv_for( 1 , 0 , 1 );

shoulder_t =0;

max_e_e = muscles1.an.pcsa*muscles1.an.m_arm_e*vars.norm_force*n_t + muscles1.tb.pcsa*muscles1.tb.m_arm_e*vars.norm_force*n_t;

max_e_f = muscles1.bb.pcsa*muscles1.bb.m_arm_e*vars.norm_force*n_t + muscles1.bs.pcsa*muscles1.bs.m_arm_e*vars.norm_force*n_t +...
    muscles1.br.pcsa*muscles1.br.m_arm_e*vars.norm_force*n_t;

max_s_e = muscles1.dp.pcsa*muscles1.dp.m_arm_s*vars.norm_force*n_t + muscles1.tb.pcsa*muscles1.tb.m_arm_s*vars.norm_force*n_t;

max_s_f = muscles1.bb.pcsa*muscles1.bb.m_arm_s*vars.norm_force*n_t + muscles1.pc.pcsa*muscles1.pc.m_arm_s*vars.norm_force*n_t +...
    muscles1.da.pcsa*muscles1.da.m_arm_s*vars.norm_force*n_t;

max_e_e
max_e_f
max_s_e
max_s_f