function [ total, mine ] = minetti( muscles , act1 , u1, vars )

for ii = 1:length(act1)
    F_max = muscles.pcsa*Fl_Fv_for(muscles.length(ii)/muscles.l0,0,1);
    a = act1(ii);
    
    v_CE_norm = muscles.v(ii)/muscles.l0;
    v_CE_max = 12;
    
    phi = (0.054+0.506*v_CE_norm+2.46*muscles.v(ii)^2)/...
        (1-1.13*v_CE_norm+12.8*(v_CE_norm^2)-1.64*(v_CE_norm)^3);
    
    Edot(ii) = a*F_max*v_CE_max*phi;
end

total =  Edot*vars.time_inc;
mine.total = total;
end

