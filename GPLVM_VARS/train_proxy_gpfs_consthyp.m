function [hyp2, datastruc,inp_ind] = train_proxy_gpfs_consthyp(input_cmg_mat,out_data,tar_data,cov_id_cell_fcst)

% this function is same as train_proxy_gpfs.m but it takes additioanl input
% of covid_mat, which gives the matrix of optimum covariance functions for
% each output for last 10 iterations. From these, it can be found out that
% which covariance is perfomring well for which output values. In this way,
% we don't need to search for optimum covariance function in training of
% GPFS model which saves a lot of time. 


[cter,lsf] = size(cov_id_cell_fcst);
[regnum,wer] = size(cov_id_cell_fcst{1,1});
for i = 1 : cter
    for j = 1 : regnum
    covid_mat(i,j) = cov_id_cell_fcst{i,1}{j,1};
    end
   
end

for i = 1 : regnum
    covid_mode(i,1) = mode(covid_mat(:,i));
end


log_input_cmg_mat = (input_cmg_mat);
clear T_vector2


T_vector2 = log_input_cmg_mat;

[u_tv,v_tv] = size(T_vector2);

well_num = 1;

% poly_vect = T_vector2;
% 
% count_rr = u_tv;
% for j = 1 : u_tv
%     for i = j : u_tv
%         count_rr = count_rr+1;
%         poly_vect(count_rr,:) = T_vector2(j,:).*T_vector2(i,:);
%     end
% end
%  T_vector2 = poly_vect;

[ererw,runs] = size(input_cmg_mat);


[lmv, gmv] = lmv_gmv_func(out_data, tar_data);

comb_inp_init{1,1} = lmv;  % LMV are the outputs of GP model
[comb_inp{1,1},inp_ind] = removeconstantrows(comb_inp_init{1,1});
% comb_inp{1,1} = mean_total_error_single ;
% comb_inp{2,1} = mean_oil_error_single ;
% comb_inp{3,1} = mean_water_error_single ;
% comb_inp{4,1} = mean_gas_error_single ;
u=1;
P1 = cell(u,1);
p1rows = cell(u,1);
P1n = cell(u,1);
p1s = cell(u,1);


for i = 1 : u
    
    
    [P1{i,1}, p1rows{i,1}] =  removeconstantrows(comb_inp{i,1});
    [P1n{i,1},p1s{i,1}] = mapminmax(P1{i,1},0,1);
end

mean_pn = mean(P1n{i,1},2);
P1n{i,1} = gsubtract(P1n{i,1},mean_pn);

global Target_cell
Target_cell = cell(1,1);
Target_cell = P1n;

T1 = cell(u,1);
t1rows = cell(u,1);
T1n = cell(u,1);
t1s = cell(u,1);


global Input_cell
Input_cell = cell(1,1);
Input_cell{1,1} = T_vector2;

global Input_cell_train
global Input_cell_val
global Input_cell_test
global Target_cell_train
global Target_cell_val
global Target_cell_test
global result_cell
global  feat_cell_train1
global f_train1
feat_cell_train1 = cell(well_num,1);
f_train1 = cell(well_num,1);

global  feat_cell_test1
global f_test1
global feat_train_2
global feat_test_2
global feat_val_2
feat_train_2 = cell(1,1);
feat_test_2 = cell(1,1);
feat_val_2 = cell(1,1);
feat_cell_test1 = cell(well_num,1);
f_test1 = cell(well_num,1);

global count
count = 0;

result_cell = cell(1000,5);

Input_cell_train = cell(u,1);
Input_cell_val = cell(u,1);
Input_cell_test = cell(u,1);

Target_cell_train = cell(u,1);
Target_cell_val = cell(u,1);
Target_cell_test = cell(u,1);

[T1{1,1}, t1rows{1,1}] =  removeconstantrows(Input_cell{1,1});
[T1n{1,1},t1s{1,1}] = mapminmax(T1{1,1},0,1);

mean_T1n = mean(T1n{1,1},2);

T1n{1,1} = gsubtract(T1n{1,1},mean_T1n);

Input_cell = T1n;

inpt_vect = (Input_cell{1,1});


for iter = 1 : 1
    
    for data_div = 1 : 1
        [inpt_train,inpt_val,inpt_test,trainInd,valInd,testInd] = dividerand(inpt_vect,0.80,0.0,0.20);
        
        mean_train = mean(inpt_train,2);
        mean_test = mean(inpt_test,2);
        
        diff_mean(data_div,1) = mean(abs(mean_train-mean_test));
        ind_store_train{data_div,1} = trainInd;
        ind_store_val{data_div,1} = valInd;
        ind_store_test{data_div,1} = testInd;
        
    end
    
    [min_val,min_loc] = min(diff_mean);
    
    trainInd = ind_store_train{min_loc,1};
    valInd = ind_store_val{min_loc,1};
    testInd = ind_store_test{min_loc,1};
    
    [inpt_train,inpt_test] = divideind(inpt_vect,trainInd,testInd);
    [T_train,T_test] = divideind(Target_cell{1,1},trainInd,testInd);
    
    x = inpt_train';
    y = T_train';
    xs = inpt_test';
    ys = T_test';
    
    [runs,num_dim] = size(y);
    
    for jj = 1 : num_dim
    [hyp2{jj,1}, mu{jj,1}, s2{jj,1}, datastruc{jj,1}] = train_gp_mult_constcov(x,y(:,jj),xs,covid_mode(jj,1));
     sdev(jj,:) = 2.*sqrt(s2{jj,1});
    mu_store(jj,:) = mu{jj,1}';
    end
    
    lb_test = mu_store - sdev;
    ub_test = mu_store + sdev;

figure(1111)
plotregression(T_test,mu_store,'GPFS Model Performance');
hold on
for jj = 1 : num_dim
scatter(T_test(jj,:),lb_test(jj,:),'r','filled');
hold on
scatter(T_test(jj,:),ub_test(jj,:),'r','filled');
hold on
end
legend('Y=X line','Regression Fit','Data','95% Credible Interval')
xlabel('Target Values')
ylabel('GP Model Output')
ax = gca;
ax.FontSize = 14;
hold off
drawnow


save next_gen_model.mat
save forward_sol.mat






end