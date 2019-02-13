% function [ shoulder, elbow, theta ,eff_mass] = looper( Data, forearm , upperarm, vars, time_step )
function [ shoulder, elbow, theta, muscles, act, u, est , tnew , energy , eff_mass ]...
    = looper2( Data, forearm , upperarm, vars, time_step,c,subj,s,t)
    %% Inverse dyanmics to joint torques
    [ shoulder , elbow , theta , eff_mass] = ...
        Calc_kine( Data, forearm, upperarm,vars.time_inc,vars.masses);

    %% Minimization of cost function by time step
    
    [ muscles , act , u , est , tnew] =...
        compute_min_cost(shoulder,elbow,theta,...
        upperarm,forearm,vars,time_step,Data,c,subj,s,t);
        
    %% Calc Energy Rates
    muscle_nums = {'an','bs','br','da','dp','pc','bb','tb'};
    for k=1:8
        [~,energy.(muscle_nums{k})] = ...
            umberger2(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}) );
    end
    
end

