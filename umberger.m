function [ total, energy ] = umberger( muscles , act , u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% muscle_nums = {'one','two','thr','fou','fiv','six'};

for ii = 1:length(act)
    F_iso = muscles.pcsa*Fl_Fv_for(muscles.length(ii),0,act(ii));
    if u(ii) > act(ii)
        A_AM = u(ii).^0.6;
    else
        A_AM = (u(ii) + act(ii)).^0.6;
    end
    
    S = 1.25;
    
    if muscles.length(ii) <= muscles.l0
        h_AM(ii) = A_AM*S*(1.28*muscles.ft+25);
    else
        h_AM(ii) = A_AM*S*(1.28*muscles.ft+25)*(.4+.6*F_iso);
    end
    
    v_CE_norm = muscles.v(ii)/muscles.l0;
    
    if u(ii)>act(ii)
        A=u(ii);
    else
        A = (u(ii)+act(ii))/2;
    end
    
    if muscles.v >= 0
        
        a_ST = 100/(12/2.5);
        a_FT = 153/12;
        
        v_CE_norm = muscles.v(ii)/muscles.l0;
        
        if muscles.length(ii) <= muscles.l0
            h_SL(ii) = -(a_ST*v_CE_norm*(1-muscles.ft)+a_FT*v_CE_norm*muscles.ft)*A^2*1.5;
        else
            h_SL(ii) = -(a_ST*v_CE_norm*(1-muscles.ft)+a_FT*v_CE_norm*muscles.ft)*A^2*1.5*F_iso;
        end
        
        w(ii) = muscles.force(ii)*muscles.v(ii);
    else
        a_L = 0.3*100/(12/2.5);
        if muscles.length(ii)<=muscles.l0
            h_SL(ii) = a_L * v_CE_norm * A * S;
        else
            h_SL(ii) = a_L * v_CE_norm * A * S * F_iso;
        end
        w(ii) = 0;
    end
    
    
end
    
    energy.h_SL = h_SL*.0025;
    energy.h_AM = h_AM*.0025;
    energy.w = w*.0025 / muscles.m;
    
    total = (energy.h_SL + energy.h_AM + energy.w)*muscles.m;
    if ~isreal(total)
        1;
    end
end

