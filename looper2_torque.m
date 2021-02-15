function [ shoulder, elbow, theta, eff_mass ]...
    = looper2_torque( Data, forearm , upperarm, vars)
    %% Inverse dyanmics to joint torques
    [ shoulder , elbow , theta , eff_mass] = ...
        Calc_kine( Data, forearm, upperarm,vars.time_inc,vars.masses);

end