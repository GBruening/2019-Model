% function [ shoulder, elbow, theta ,eff_mass] = looper( Data, forearm , upperarm, vars, time_step )
function [ shoulder, elbow, theta, muscles, act, u, est , tnew ,...
    umber, bhar, uch, lich, marg, mine, houd, eff_mass ]...
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
    % Calc houdijk Energy Rates
        [~,houd.(muscle_nums{k})] = ...
            houdijk(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars );
        
    % Calc Minetti Energy Rates
        [~,mine.(muscle_nums{k})] = ...
            minetti(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars  );
        
    % Calc Umberger Energy Rates
        [~,umber.(muscle_nums{k})] = ...
            umberger2(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}) );
        
    % Calc Bhargava Energy Rates
        [~,bhar.(muscle_nums{k})] = ...
            bhargava(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars );
        
    % Calc uchida Energy Rates
        [~,uch.(muscle_nums{k})] = ...
            uchida(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars );
        
    % Calc lichtwark Energy Rates
        [~,lich.(muscle_nums{k})] = ...
            lichtwark(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars );
        
    % Calc Margaria Energy Rates
        [~,marg.(muscle_nums{k})] = ...
            margaria(muscles.(muscle_nums{k}) , act.(muscle_nums{k}) , u.(muscle_nums{k}), vars );
        
    end
end

