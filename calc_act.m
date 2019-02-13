function [ act ] = calc_act( muscles , ii )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    muscle_nums = {'one','two','thr','fou','fiv','six'};
    
    
    
    if ii~=length(muscles.one.m_arm)
        for k=1:length(muscle_nums)
            norm_length = muscles.(muscle_nums{k}).length(ii)/muscles.(muscle_names{k}).l0;
            vel = muscles.(muscles_nums{k}).v(ii);
            act.(muscle_nums{k})(ii) = Fl_Fv(norm_length,vel,
        end
    else

end

