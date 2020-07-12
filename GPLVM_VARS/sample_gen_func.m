function [input_cmg_mat] = sample_gen_func(input_fwd_sol,Target_pred_result,Target_fwd_pred_result)

% this function joins columns of three matrices 

input_cmg_mat = input_fwd_sol;
[u_t,v_t] = size(input_cmg_mat);

[u_ss,v_ss] = size(Target_pred_result);
[u_ss1,v_ss1] = size(Target_fwd_pred_result);
if u_ss == 0 && u_ss1 == 0
elseif u_ss ~= 0    
input_cmg_mat(:,v_t+1 : v_t+v_ss ) = Target_pred_result(:,1 : v_ss);
elseif u_ss1 ~= 0
    [u_t,v_t] = size(input_cmg_mat);
input_cmg_mat(:,v_t+1 : v_t+v_ss1) = Target_fwd_pred_result(:,1 : v_ss1);
end


end