function [lb,ub,lb_1,ub_1] = calc_bound_func(perc_inp_mat,num_out_cell,error_mat,input_cmg_mat, mean_input_mat,count_er,change_par,lower_thers,upper_thers,max_area_final,max_area_ind_final)

%% This function calculates the updated lower and upper bounds of each input parameter
%outputs 
% lb = lower bound
% ub = upper bound
% lb_1 = not important, please ignore
% ub_1 = not important, please ignore
%%

[u_t,v_t] = size(input_cmg_mat);
[j_t,j_m] = size(error_mat);
% [u_ch,v_ch] = size(change_par);

norm_error_mat = abs(error_mat./max(abs(error_mat)));
wt_error_mat = norm_error_mat./ sum(norm_error_mat);

for ii = 1 : u_t
    check_sign = 0;
    err_wt = 0;
    for jj = 1 : j_t
        check_sign_cell{ii,1}(jj,1) = num_out_cell{ii,1}(jj,1)/num_out_cell{ii,2}(jj,1);
        if check_sign_cell{ii,1}(jj,1)<0
             check_sign = check_sign+1;
             err_wt = err_wt + wt_error_mat(jj,1);
        end
        
        
    end
    err_wt_mat(ii,1) = err_wt;
    
    
    
%     if check_sign==j_t
%         sign_mat(ii,1) = check_sign;
%     else
%         sign_mat(ii,1) = 0;
%     end
%      if check_sign==j_t-1
%         sign_mat(ii,1) = check_sign;
%      end
         sign_mat(ii,1) = check_sign;
end
        
        
check_err = 0;
for i = 1 : u_t
    check_num1 = 0;
    check_num2 = 0;
    for j = 1 : j_t
        if abs(error_mat(j,1))<1
            check_err = 1;
        end
            
        
        if error_mat(j,1)>0
            if num_out_cell{i,1}(j,1)<0
                bound_cell{i,1}(j,1) = (-abs((perc_inp_mat(i,1)/num_out_cell{i,1}(j,1))*error_mat(j,1))) * max_area_ind_final(i,j);
                %                 check_num1 = check_num1+1;
            elseif num_out_cell{i,2}(j,1)<0
                bound_cell{i,1}(j,2) = (abs((perc_inp_mat(i,2)/num_out_cell{i,2}(j,1))*error_mat(j,1))) * max_area_ind_final(i,j);
                %                 check_num2 = check_num2+1;
            else
                check_num1 = check_num1+1;
                check_num2 = check_num2+1;
                
            end
            
        end
        if error_mat(j,1)<=0
            if num_out_cell{i,1}(j,1)>0
                bound_cell{i,1}(j,1) = (-abs(perc_inp_mat(i,1)/(num_out_cell{i,1}(j,1))*error_mat(j,1)))* max_area_ind_final(i,j);
                %                  check_num1 = check_num1+1;
            elseif num_out_cell{i,2}(j,1)>0
                bound_cell{i,1}(j,2) = (abs((perc_inp_mat(i,2)/num_out_cell{i,2}(j,1))*error_mat(j,1)))* max_area_ind_final(i,j);
                %                 check_num2 = check_num2+1;
            else
                check_num1 = check_num1+1;
                check_num2 = check_num2+1;
                
            end
        end
        
        if check_num1==j_t
            bound_cell{i,1}(j,1) = 0;
        end
        if check_num2==j_t
            bound_cell{i,1}(j,2) = 0;
        end
        if num_out_cell{i,1}(j,1)==0
            bound_cell{i,1}(j,1) = 0;
            
        end
        if num_out_cell{i,2}(j,1)== 0
            bound_cell{i,1}(j,2) = 0;
        end
    end
    perc_bd_mat(i,1) = min(min(bound_cell{i,1}));
    
      if abs(perc_bd_mat(i,1)) <= 0.05
        perc_bd_mat(i,1) = -0.20;
      end
    
    if  perc_bd_mat(i,1)<= -1
        if (max_area_final(i,1)>=0.9)
        perc_bd_mat(i,1) = -max_area_final(i,1);
        else
           perc_bd_mat(i,1) = -max_area_final(i,1);
        end
    end
    perc_bd_mat(i,2) = max(max(bound_cell{i,1}));
    
    if perc_bd_mat(i,2) <= 0.05
        perc_bd_mat(i,2) = 0.20;
    end
    
    if  perc_bd_mat(i,2)>= 1
          if (max_area_final(i,1)>=0.9)
        perc_bd_mat(i,2) = max_area_final(i,1);
        else
           perc_bd_mat(i,2) = max_area_final(i,1);
        end
    end
    
%     check_ch = 0;
%     for jj = 1 : u_ch
%         if i == change_par(jj,1)
%             check_ch = 1;
%             break
%         end
%     end
%     
%     if check_ch == 0
%          perc_bd_mat(i,1) = 0;
%           perc_bd_mat(i,2) = 0;
%     end
    
%     if   (sign_mat(i,1)==0)
%         perc_bd_mat(i,1)= 0;
%         perc_bd_mat(i,2)= 0;
%     end
%     
%     if i<=3
%         if (check_num1>0) || (check_num2>0) || (sign_mat(i,1)<j_t)
%             perc_bd_mat(i,1)= 0;
%             perc_bd_mat(i,2)= 0;
%         end
%     end
%     if i==u_t
%         if (check_num1>1) || (check_num2>1)|| (sign_mat(i,1)<j_t-1)
%             perc_bd_mat(i,1)= 0;
%             perc_bd_mat(i,2)= 0;
%         end
%     end
 
  
   
   
    
    
%     if (check_err==1) && (count_er>=3)
%         perc_bd_mat(u_t,1)= -0.05;
%          perc_bd_mat(u_t,2)= 0.05;
%     end
%  


end

for i = 1 : u_t
   sig = perc_bd_mat(i,2)/ perc_bd_mat(i,1);
    
    if sig < 0

       tot_perc(i,1) = perc_bd_mat(i,2) - perc_bd_mat(i,1);
    else
        tot_perc(i,1) = max(max(abs(perc_bd_mat(:,2) - perc_bd_mat(:,1))));
    end
end
norm_tot_perc = zeros(u_t,1);
for i = 1 : u_t
    if tot_perc(i,1)~= 0;
norm_tot_perc(i,1) = tot_perc(i,1)./max(abs(tot_perc));
    else
        norm_tot_perc(i,1) = 0;
    end
end

for i = 1 : u_t
    if norm_tot_perc(i,1) ~= 0

inv_tot_perc(i,1)= 1/ norm_tot_perc(i,1);
    else
        inv_tot_perc(i,1)= 0;
    end
end


norm_tot_perc = max_area_final;



% norm_tot_perc = norm_tot_perc./ max(norm_tot_perc);

norm_sign_mat = sign_mat./max(sign_mat);
count_ch = 0;
for i = 1 : u_t
%     if sign_mat(i,1) < int64(max(sign_mat)/2)
%     if norm_tot_perc(i,1) <= 0.5
%         norm_tot_perc(i,1) = 0.1;
%     end
%     end
    fact(i,1) = 1.0;
%     if fact(i,1)>1
%         fact(i,1) =1;
%     end
    lb(i,1) = input_cmg_mat(i,1) - abs(perc_bd_mat(i,1)*fact(i,1)*(mean_input_mat(i,1)));
    ub(i,1) = input_cmg_mat(i,1) + abs(perc_bd_mat(i,2)*fact(i,1)*(mean_input_mat(i,1)));
%     lb(i,1) = input_cmg_mat(i,1) - abs(perc_bd_mat(i,1)*(mean_input_mat(i,1)));
%     ub(i,1) = input_cmg_mat(i,1) + abs(perc_bd_mat(i,2)*(mean_input_mat(i,1)));

for jj = 1 : u_t
%     if i == change_par(jj,1)
    count_ch = count_ch+1;
    lb_1(count_ch,1) = lower_thers(i,1);
    ub_1(count_ch,1) = upper_thers(i,1);
%     end
end


end


end