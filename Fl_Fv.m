function [ a ] = Fl_Fv( l , v, F )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%     Fp1 = 68.35*0.0495*log(exp((l-1.445)/(0.0495))+1);
    Fp2 = -0.02*exp(13.8-18.7*l);
    
    if v <= 0
        Fv = (-5.72-v)/(-5.72+(1.38+2.09*l)*v);
    else
        Fv = (0.62-(-3.12+4.21*l-2.67*(l^2))*v)/(0.62+v);
    end

    Fl = exp(-(abs((l^1.93-1)/1.03)^1.87));

    Nf = 2.11 + 4.16*((1/l)-1);
    
    syms x
%     a = solve( F == (1-exp(-((x/(0.56*Nf))^Nf)))*(Fl*Fv+Fp2)+Fp1,x);
    a = solve( F == (1-exp(-((x/(0.56*Nf))^Nf)))*(Fl*Fv+Fp2),x);
    a = double(a);
    
end

