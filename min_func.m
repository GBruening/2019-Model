function [output] = min_func(x,muscles,act,vars,ii,Aeq,beq)

if  isa(act,'struct')
else
    clear act
end

muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};

time_step = .0025;
output_stress = sum(x);
output_stress2 = sum(x.^2);

switch vars.minparam
    case 'stress'
        output = output_stress;
        return
    case 'stress2'
        output = output_stress2;
        return
end

muscles.an.force(ii,1) = x(1)*muscles.an.pcsa;
muscles.bs.force(ii,1) = x(2)*muscles.bs.pcsa;
muscles.br.force(ii,1) = x(3)*muscles.br.pcsa;
muscles.da.force(ii,1) = x(4)*muscles.da.pcsa;
muscles.dp.force(ii,1) = x(5)*muscles.dp.pcsa;
muscles.pc.force(ii,1) = x(6)*muscles.pc.pcsa;
muscles.bb.force(ii,1) = x(7)*muscles.bb.pcsa;
muscles.tb.force(ii,1) = x(8)*muscles.tb.pcsa;

for muscle = muscle_nums
    F_max = muscles.(muscle{1}).pcsa*Fl_Fv_for(muscles.(muscle{1}).length(ii)/muscles.(muscle{1}).l0,0,1);
    if muscles.(muscle{1}).force(ii,1)/F_max>10
        1;
    end
end

min_force = 1E-10;
output_force = 0;
output_force2 = 0;
for k=1:length(muscle_nums)
    if muscles.(muscle_nums{k}).force(ii) < 0
        muscles.(muscle_nums{k}).force(ii) = 0;
    end
    output_force = output_force + muscles.(muscle_nums{k}).force(ii,1);
    output_force2 = output_force2 + muscles.(muscle_nums{k}).force(ii,1)^2;
end

switch vars.minparam
    case 'force'
        output = output_force;
        return
    case 'force2'
        output = output_force2;
        return
end

norm_force = vars.norm_force; %31.8E4
for k=1:length(muscle_nums)
    norm_length = muscles.(muscle_nums{k}).length(ii)/muscles.(muscle_nums{k}).l0;
    vel         = muscles.(muscle_nums{k}).v(ii);
    n_f         = muscles.(muscle_nums{k}).force(ii)/(norm_force*muscles.(muscle_nums{k}).pcsa);    
    act.(muscle_nums{k})(ii) = Fl_Fv_inv(norm_length,vel,n_f);
end

output_act = 0;
output_act2 = 0;
for k=1:length(muscle_nums)
    output_act  = output_act  + act.(muscle_nums{k})(ii);
    output_act2 = output_act2 + act.(muscle_nums{k})(ii).^2;
end

switch vars.minparam
    case 'act'
        output = output_act;
        return
    case 'act2'
        output = output_act2;
        return
end

t_act = 0.05;%.050;
t_deact = 0.066;%.066;

sig = 0;
eps = 0;

for k=1:length(muscle_nums)
    if ii==1

        adot = act.(muscle_nums{k})(1)/time_step;
        
        u1 = (adot*t_deact + 0)./(1+sig*eps-adot*(t_act-t_deact));
        u2 = (adot*t_deact + 0)./(1+sig*eps);
        
        if u1>act.(muscle_nums{k})(1)
            drive(k) = u1;
        else
            drive(k) = u2;
        end
    else
        adot = (act.(muscle_nums{k})(ii)-act.(muscle_nums{k})(ii-1))/time_step;
        u1   = (adot*t_deact + act.(muscle_nums{k})(ii-1))./(1+sig*eps-adot*(t_act-t_deact));
        u2   = (adot*t_deact + act.(muscle_nums{k})(ii-1))./(1+sig*eps);
        
        if u1>act.(muscle_nums{k})(ii)
            drive(k) = u1;
        else
            drive(k) = u2;
        end
    end
    drive(k) = drive(k);
    scale_drive(k) = drive(k)*muscles.(muscle_nums{k}).m;
end

output_drive = sum(scale_drive);
output_drive2 = sum(scale_drive.^2);

switch vars.minparam
    case 'drive'
        output = output_drive;
        return
    case 'drive2'
        output = output_drive2;
        return
end


for k=1:length(muscle_nums)
    [energy(k),total] = umberger2(muscles.(muscle_nums{k}),act.(muscle_nums{k})(ii),drive(k));
end
    
if sum(energy) == inf || isnan(sum(energy))
    1;
end
    
switch vars.minparam
    case 'umberger'
        output = sum(energy);
        return
end

end
