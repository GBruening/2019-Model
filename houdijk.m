function [ total, houd ] = houdijk( muscles , act1 , u1, vars )

hbar_af = 52.5; %W/kg
hbar_as = 10.98; %W/kg

hbar_mf = 97.5; %W/kg
hbar_ms = 13.42; %W/kg

k1_f = 12;
k1_s = 6;

k2_f = 14;
k2_s = 8;


for ii = 1:length(act1)
    F_iso = muscles.pcsa*Fl_Fv_for(muscles.length(ii)/muscles.l0,0,act1(ii));
    
    F_max = muscles.pcsa*Fl_Fv_for(muscles.length(ii)/muscles.l0,0,1);
    
    F = muscles.force(ii)/F_max;
    
    hbar_slf = 0.28*F_max;
    hbar_sls = 0.16*F_max;
    
    vmax_f = k1_f+k2_f*act1(ii);
    vmax_s = k1_s+k2_s*act1(ii);
    
    ha_f(ii) = muscles.ft*muscles.m*hbar_af*u1(ii)*...
        (1-exp(-0.25-(18.2/(u1(ii)*vmax_f))))/...
        (1-exp(-0.25-(18.2/vmax_f)));
    ha_s(ii) = (1-muscles.ft)*muscles.m*hbar_as*u1(ii)*...
        (1-exp(-0.25-(18.2/(u1(ii)*vmax_s))))/...
        (1-exp(-0.25-(18.2/vmax_s)));
    
    hm_f(ii) = muscles.ft*muscles.m*(hbar_af+hbar_mf)*act1(ii)*(F-(hbar_af)/(hbar_af+hbar_mf));
    hm_s(ii) = (1-muscles.ft)*muscles.m*(hbar_as+hbar_ms)*act1(ii)*(F-(hbar_as)/(hbar_as+hbar_ms));

    hsl_f(ii) = muscles.ft*hbar_slf*act1(ii)*F*muscles.v(ii);
    hsl_s(ii) = (1-muscles.ft)*hbar_slf*act1(ii)*F*muscles.v(ii);
    
%     G = muscles.1
    
    v_CE_norm = muscles.v(ii)/muscles.l0;
    if v_CE_norm>=0
        Wdot(ii) = (muscles.force(ii)*muscles.v(ii));
    else
        Wdot(ii) = 0;
    end
end

houd.h_SL = (abs(hsl_f)+abs(hsl_s)) * vars.time_inc;
houd.h_AM = ((abs(ha_f)+abs(ha_s))+(abs(hm_f)+abs(hm_s))) * vars.time_inc * muscles.m;
houd.w = abs(Wdot) * vars.time_inc;

total = (houd.h_SL + houd.h_AM + houd.w);
houd.total = total;
end

