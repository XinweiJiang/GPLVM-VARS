close all
% clc
clear all

%% Before running this code, please give the path where GPML, GPLVM and Variogram matlab packages are installed below

path_gpml = '..\gpml-matlab-v4.0-2016-10-19';                     % please insert here path of GPML Toolbox
path_vario = '..\variogram';                              % please insert here path of Variogram Toolbox
path_gplvm = '..\gplvm';                              % please insert here path of GPLVM Toolbox

% bUseGPLVM4GPIS: 0/1 (train_proxy_gpis with GPR(0) or GPLVM(1))
% 0: using original GP_VARS and GPIS in train_proxy_gpis
% 1: using the proposed GPLVM-VARS and GPLVM based GPLVMIS in train_proxy_gpis
bUseGPLVM4GPIS = 1;

% add folders of MTGP and GPML Toolbox
if ~isunix  % windows system
    if bUseGPLVM4GPIS,   addpath(genpath(path_gplvm));   end;
    addpath(genpath(path_gpml));
    addpath(genpath(path_vario));
else        % linux system
    addpath(genpath('../../'));
    if bUseGPLVM4GPIS,   addpath(genpath(path_gplvm));   end;
    addpath(genpath(path_gpml));
    addpath(genpath(path_vario));
end
%%

%% Information from the user on which function to test on   

fun = 2;   % 1 for running paper example problem. 2 for running Rastrigin function problem



lower_thers = xlsread('input_data.xlsx',['lower_limit_case' num2str(fun)]);
upper_thers = xlsread('input_data.xlsx',['upper_limit_case' num2str(fun)]);

runs = 50;  % initial number Latin Hypercube Samples

%% Prepare training data from initial input and output data

[input_cmg_mat] = lhs_func_rev1(lower_thers,upper_thers,[],[],runs);  % Initial Latin hypercube samples

for i = 1 : runs
    if fun ==1
        out_data(:,i) = example_prob(input_cmg_mat(:,i));  % Obtain output corresponding to selected inputs
    else
        out_data(:,i) = Rastrigin(input_cmg_mat(:,i));
    end
    
end

tar_data = xlsread('input_data.xlsx',['target_data_case' num2str(fun)]);    % This is the target observed value 

%%

comb_inp_error =100000;    % This is same as GMV

 lb = lower_thers;   % lower limit of input parameters
 ub = upper_thers;   % upper limit of input parameters

 orig_input_cmg_mat = input_cmg_mat;
 orig_out_data = out_data;
 orig_tar_data = tar_data;
 
 %% Train inverse gp model from training data 
[hyp_inv,datastruc_inv,cov_id_cell_inv{1,1}] = train_proxy_gpis(input_cmg_mat,out_data,tar_data, bUseGPLVM4GPIS);
    
 [Target_pred_result] = func_proposal_sample_inverse_gp();  % Obtain first proposal solution
    
str{1,1} = 'next_gen_model.mat';
inv_ann_result = Target_pred_result; % storing the first proposal sample
input_par = inv_ann_result;


save upper_half_const.mat    % save all the data so far

b_check = 0;
count_rec = 1;
uorig_case = 1;
global ep
ep = 0.01;

global count_er
count_er=0;      % count number of iterations for GP-VARS

global tol_bayes
tol_bayes = 0.01;

T_vect_rev = input_cmg_mat;    % T_vect_rev will contain all training data that be acquired during GP-VARS iterations

next_out_data = orig_out_data;  % next_out_data will contain all training output data that be acquired during GP-VARS iterations

tttime = cputime;

% while(comb_inp_error>0.001 &&  count_er<100)
while(comb_inp_error>0.001)
    count_er = count_er + 1;
  
    str{1,1} = 'next_gen_model.mat';
    
    if count_er == 1
        Target_pred_result = inv_ann_result;
    else
        
        Target_pred_result = pred_result_best(:,count_er-1)
        input_par = Target_pred_result;
        
    end
    
    clear orig_input1   
    clear input_cmg_mat
    
    input_cmg_mat = Target_pred_result;
    
    
 
    saved_result_cell{count_er,1} = input_cmg_mat;
  
    
    [u_t,v_t] = size(input_cmg_mat);
    clear result_out
 
    for i = 1 : v_t
        if fun ==1
            result_out(:,i) = example_prob(input_cmg_mat(:,i));  % Obtain output corresponding to selected inputs
        else
            result_out(:,i) = Rastrigin(input_cmg_mat(:,i));
        end
    
end
    
    [lmv, gmv] = lmv_gmv_func(result_out, tar_data);
    
           comb_inp_error = gmv;
           [error_mat] = lmv;

      saved_out_data{count_er,1} = result_out;
    
    mod_error_mat(count_er,1) = count_er;
    mod_error_mat(count_er,2) = comb_inp_error;
    
    
    
    %%    Train forward GP Model (GPFS)
    if count_er == 1
        
        [hyp{i,1}, datastruc{i,1},inp_ind_fwd_mat{i,1},cov_id_cell_loop{i,1}] = train_proxy_gpfs(orig_input_cmg_mat,orig_out_data,orig_tar_data);
        
    else
        if count_er <= 10
            [hyp{i,1}, datastruc{i,1},inp_ind_fwd_mat{i,1},cov_id_cell_loop{i,1}] = train_proxy_gpfs(T_vect_rev,next_out_data,orig_tar_data);
        else
            [hyp{i,1}, datastruc{i,1},inp_ind_fwd_mat{i,1}] = train_proxy_gpfs_consthyp(T_vect_rev,next_out_data,orig_tar_data,cov_id_cell_fwd);
        end
        
        end
    %% store the id of covariance function that was the most optimum in training of GPFS
    inp_ind_fwd = inp_ind_fwd_mat{1,1};
    hyp_fwd{count_er,1} = hyp{1,1};
    datastruc_fwd{count_er,1} = datastruc{1,1};
    cov_id_cell_fwd{count_er,1} = cov_id_cell_loop{1,1};
    
     
    % check error for forward ann
    
    %%
    % Run VARS based sensitivity analysis
        [max_area_final,max_area_ind_final] = func_vars_sensitivity1();

    %% Calculate updated lower and upper bounds
    [perc_inp_mat,num_out_cell,mean_input_mat,cell_sens_resp] = sens_ana_func_gp(lower_thers,upper_thers);

   [lb,ub,lb_1,ub_1] = calc_bound_func(perc_inp_mat,num_out_cell,error_mat,input_cmg_mat, mean_input_mat,count_er,[],lower_thers,upper_thers,max_area_final,max_area_ind_final);

    [u_t,v_t] = size(lb);
    
    if count_er>10
        uorig_case = uorig;
    end
    % Check if updated lower and upper bounds are outside original ranges
    % of lower and upper limits
    
    for ii = 1 : u_t
        if lb(ii,1)<lower_thers(ii,1)
            lb(ii,1)=lower_thers(ii,1);
        end
        if lb(ii,1)>upper_thers(ii,1)
            lb(ii,1)=upper_thers(ii,1);
        end
        if ub(ii,1)<lower_thers(ii,1)
            ub(ii,1)=lower_thers(ii,1);
        end
        if ub(ii,1)>upper_thers(ii,1)
            ub(ii,1)=upper_thers(ii,1);
        end
    end
    
    lb_cell{count_er,1} = lb;   % Store lb and ub for different iterations
    ub_cell{count_er,1} = ub;
    orig_input1 = orig_input_cmg_mat;
    if count_er == 1
        [u_t,v_t] = size(orig_input1);
    else
        [u_t,v_t] = size(prev_input_cmg_mat_store);
    end
    
    
        b_check=1;
        count_rec = count_er;
        uorig_case = u_t;

    
    
    if (count_er>4) && (v_t>30) && (mod_error_mat(count_er,2)>0.05)
        if  ((mod_error_mat(count_er,2) + mod_error_mat(count_er-1,2)> mod_error_mat(count_er-2,2) + mod_error_mat(count_er-3,2))) ||  (abs(mod_error_mat(count_er,2)-mod_error_mat(count_er-2,2))<0.01)
            
            b_check=1;
            count_rec = count_er;
            uorig_case = uorig-1;
        else
            b_check=0;
        end
    else
        b_check=0;
    end
    
    if  (count_er>2) && (abs(mod_error_mat(count_er,2)-mod_error_mat(count_er-2,2))<=0.005)
        b_check=1;
        count_rec = count_er;
        uorig_case = uorig;
    end
    

    if  b_check==0
        lb = lb_cell{count_rec,1};
        ub = ub_cell{count_rec,1};
    end
    
    if count_er>2
        [u_t,v_t] = size(lb);
        for i = 1 : u_t

            if (ub(i,1)-lb(i,1)) > (ub_cell{count_er-1,1}(i,1)-lb_cell{count_er-1,1}(i,1))
                diff_b = (ub_cell{count_er-1,1}(i,1)-lb_cell{count_er-1,1}(i,1));
                ub(i,1) = ub(i,1) - abs(diff_b/4);
                lb(i,1) = lb(i,1) + abs(diff_b/4);
            end
            
            if (ub(i,1)-lb(i,1))/lb(i,1) < 0.01
                diff_b = (ub_cell{count_er-1,1}(i,1)-lb_cell{count_er-1,1}(i,1));
                ub(i,1) = ub(i,1) + abs(diff_b/4);
                lb(i,1) = lb(i,1) - abs(diff_b/4);
            end
        end
    end

    lb_cell{count_er,1} = lb;
    ub_cell{count_er,1} = ub;
    
        save mid_code.mat
%%  Generate 10 samples via Bayesian Optimization

    [input_vect_store] = func_proposal_sample_bayes_opt(str,lb,ub);

    Target_pred_input = input_cmg_mat;
        
    

    tol = 0.02;
    count_next = 0;
    repeat_qty = 0;
    if (count_er >=1) || (check_flg > 0)
        
        if (count_er>1)
            [u_n,v_n] = size(T_vect_rev);
            [u_n1,v_n1] = size(prev_input_cmg_mat_store);
            for ii = 1 : v_n1
                diff_mat = gsubtract(T_vect_rev,prev_input_cmg_mat_store(:,ii));
                [uin,vin] = size(diff_mat);
                for coln = 1 : vin
                norm_diff_mat(1,coln) = norm(diff_mat(:,coln));
                end
                
                if min(norm_diff_mat) > tol
                     count_next = count_next + 1;
                
                T_vect_rev(:,count_next+v_n) = prev_input_cmg_mat_store(:,ii);
                next_out_data(:,count_next+v_n) = prev_out_data_store(:,ii);
               
                else
                    repeat_qty = repeat_qty + 1;
                end
            end
            
        end
        
        lb1=lb;
        ub1=ub;
        [uorig,vorig] = size(T_vect_rev);
        found_inp = 0;
        for j = 1 : vorig
            check_lt = 0;
            for i = 1 : uorig
                if (T_vect_rev(i,j)>=lb1(i,1)) && (T_vect_rev(i,j)<=ub1(i,1))
                    check_lt = check_lt+1;
                end
            end
            if check_lt >= uorig_case
                found_inp = found_inp + 1;
            end
        end
        
        clear more_out_data
        
  %% Obtain model output for 10 proposal samples obtained via Bayesian Optimization
                
            [input_cmg_mat] = sample_gen_func(input_vect_store,[],[]);
            [wer,v_t] = size(input_cmg_mat);
 
            for i = 1 : v_t
                if fun ==1 
                more_out_data(:,i) = example_prob(input_cmg_mat(:,i));
                else
                   more_out_data(:,i) = Rastrigin(input_cmg_mat(:,i));
                end
                
            end
            
            
 %% Add the proposal samples and their true outputs to previous set of training data
            prev_input_cmg_mat_store = input_cmg_mat;
            prev_out_data_store =  more_out_data;
            
            save saved_first_round.mat

        
        [u_inp,v_inp] = size(input_vect_store);
        
        input_cmg_mat_recent = input_cmg_mat;
        out_data_inv = more_out_data;
        
        
        
        [u_t,v_t] = size(input_cmg_mat);
        [uorig,vorig] = size(T_vect_rev);
        found_inp = 0;
        
         lb1 = lower_thers;
         ub1 = upper_thers;
         
    tol = 0.02;
    count_next = 0;
    repeat_qty = 0;
        for j = 1 : vorig
            check_lt = 0;
            for i = 1 : uorig
                if (T_vect_rev(i,j)>=lb1(i,1)) && (T_vect_rev(i,j)<=ub1(i,1))
                    check_lt = check_lt+1;
                end
            end
            if count_er== 1
                if check_lt >= 1
                    found_inp = found_inp + 1;
                    input_cmg_mat(:,v_t+found_inp) = T_vect_rev(:,j);
                    out_data_inv(:,v_t+found_inp) = orig_out_data(:,j);
       
                end
            else
                if (check_lt >= uorig_case)
                    found_inp = found_inp + 1;
                end
                if found_inp>0
                    diff_mat = gsubtract(input_cmg_mat,T_vect_rev(:,j));
                [uin,vin] = size(diff_mat);
                for coln = 1 : vin
                norm_diff_mat(1,coln) = norm(diff_mat(:,coln));
                end
                
                if min(norm_diff_mat) > tol
                    count_next = count_next + 1;

                    input_cmg_mat(:,v_t+count_next) = T_vect_rev(:,j);
                    out_data_inv(:,v_t+count_next) = next_out_data(:,j);
        
                end
                end
            end
            
        end
    end
%%    Train inverse GP model (GPIS)
    
   if count_er <50
        [hyp_inv,datastruc_inv,cov_id_cell_inv{count_er,1}] = train_proxy_gpis(input_cmg_mat,out_data_inv,tar_data, bUseGPLVM4GPIS);
    else
        [hyp_inv,datastruc_inv] = train_proxy_gpis_consthyp(input_cmg_mat,out_data_inv,tar_data,cov_id_cell_inv, bUseGPLVM4GPIS);
   end
    %% Obtain 1 proposal sample via inverse GP solution (GPIS)
    
    [Target_pred_result_new] = func_proposal_sample_inverse_gp();
    Target_pred_result_new = real(Target_pred_result_new);
    [uii,vii] = size(Target_pred_result_new);
   
    % check if proposed sample via inverse GP is outside the original
    % lower and upper limits
    
    for iii = 1 : uii
        if Target_pred_result_new(iii,1) < lower_thers(iii,1)
            Target_pred_result_new(iii,1) = lower_thers(iii,1);
        end
        if Target_pred_result_new(iii,1) > upper_thers(iii,1)
            Target_pred_result_new(iii,1) = upper_thers(iii,1);
        end
    end

    %% Obtain true output corresponding to 1 proposal sample via inverse GP model (GPIS)
    for i = 1 : 1
        if fun ==1 
        more_out_data(:,v_t+i) = example_prob(Target_pred_result_new);
        else
            more_out_data(:,v_t+i) = Rastrigin(Target_pred_result_new);
        end
    end
        
    
    input_cmg_mat = sample_gen_func(input_cmg_mat_recent,Target_pred_result_new,[]);
   %% Calculate GMV for the 11 proposal samples and save all the results in the saved_result_cell variable
    [lmv, gmv] = lmv_gmv_func(more_out_data,tar_data);
    
    comb_inp_error_2 = gmv;
    comb_inp_error = min(gmv);
        
    [erwer,v_inp] = size(comb_inp_error_2);
    
    for j = 1 : v_inp
        mod_error_mat(count_er,j+2) = comb_inp_error_2(1,j);
        saved_result_cell{count_er,j+1} = input_cmg_mat(:,j);
        saved_out_data{count_er,j+1} = more_out_data(:,j);
        
    end
    
    % Best result
    [min_comb_err,min_comb_err_loc] = min(comb_inp_error_2);
    if min_comb_err > comb_inp_error
        pred_result_best(:,count_er) = Target_pred_result;
    else
        pred_result_best(:,count_er) = input_cmg_mat(:,min_comb_err_loc);
    end
    %%

    prev_input_cmg_mat_store = input_cmg_mat;
    prev_out_data_store = more_out_data;  
    
   
    save saved_next_gen_model.mat
    
    %% view the performance of all proposal samples
    
    figure(501)
    for j = 1 : v_inp
        if j<=5
            scatter(mod_error_mat(:,1),mod_error_mat(:,j+2),'k','*')
        elseif j>5 && j<=10
            scatter(mod_error_mat(:,1),mod_error_mat(:,j+2),'b','+')
        else
            scatter(mod_error_mat(:,1),mod_error_mat(:,j+2),'m','filled')
            
        end
        hold on
    end
    grid on
    xlabel('Iterations')
    ylabel('Global Misfit Value')
%     ylim([0 5])
    drawnow;
    hold off
    
     %% view results of updated lower and upper bounds
    
    
    % count_er =6;
    for i = 1 : count_er
        ct_mat(i,1) = i;
    end
    
    [u_t,v_t] = size(lb_cell{count_er,1});
    for i = 1 : u_t
        for j = 1 : count_er
            figure(i+10)
            scatter(ct_mat(j,1),lb_cell{j,1}(i,1),'b','filled')
            hold on
            scatter(ct_mat(j,1),ub_cell{j,1}(i,1),'r','filled')
            xlabel('Iterations')
            ylabel(['Bounds of Parameter' num2str(i)])
            drawnow;
        end
    end
    
    
end

e = cputime-tttime
fprintf('\n\n************   Exit while loop with comb_inp_error: %f and count_er: %d and bUseGPLVM4GPIS: %d (train_proxy_gpis with GPR(0) or GPLVM(1))    *************\n\n', comb_inp_error, count_er, bUseGPLVM4GPIS);

% load next_gen_model.mat
% load forward_sol.mat
% 
% figure(11111)
% plotregression(T_test,mu_store,'GPFS Model Performance');
% hold on
% for jj = 1 : num_dim
% scatter(T_test(jj,:),lb_test(jj,:),'r','filled');
% hold on
% scatter(T_test(jj,:),ub_test(jj,:),'r','filled');
% hold on
% end
% legend('Y=X line','Regression Fit','Data','95% Credible Interval')
% xlabel('Target Values')
% ylabel('GP Model Output')
% ax = gca;
% ax.FontSize = 14;
% hold off
% drawnow
% print(['E:\dv_GPFS_' num2str(nnnn) '.eps'],'-depsc2','-tiff', '-loose', '-r600' );
% 
% 


% best result
proposal_sample_gmv = mod_error_mat(:,2:end);
[min_val_mat,min_row] = min(proposal_sample_gmv);
[min_val,min_col] = min(min_val_mat);
best_proposal_sample = saved_result_cell{min_row,min_col};
actual_answer = xlsread('input_data.xlsx',['final_answer_case' num2str(fun)]);


