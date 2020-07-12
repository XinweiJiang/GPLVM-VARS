function [input_cmg_mat] = lhs_func_rev1(lb,ub,Target_pred_result,Target_fwd_pred_result,runs)
% this function is used to generate latin hypercube samples

lb = transpose(lb);
ub = transpose(ub);
[bb,num_par] = size(lb);
xn = lhsdesign(runs,num_par);
sample_lhs = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb)));
sample = transpose(sample_lhs);

input_cmg_mat = sample;
[u_t,v_t] = size(input_cmg_mat);

[u_ss,v_ss] = size(Target_pred_result);
[u_ss1,v_ss1] = size(Target_fwd_pred_result);
if u_ss == 0 || u_ss1 == 0
else
input_cmg_mat(:,v_t+1) = Target_pred_result(:,1);
input_cmg_mat(:,v_t+2) = Target_fwd_pred_result(:,1);
end


end