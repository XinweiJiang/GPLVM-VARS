function [hyp2,datastruc,cov_id] = train_proxy_gpis(input_cmg_mat,out_data,tar_data, bUseGPLVM4GPIS)

% This function is used to train inverse GP model. (GPIS)
% this code is almost same as the code for training GPFS model. The only
% difference is that outputs of GPFS model become inputs of GPIS model.
% inputs
% input_cmg_data = all training input data. This code standardizes this
% data later on.
% out_data = output data corresponding to the training data
% tar_data = target data or observed values

% outputs
% the main output of this function is the inverse_model.mat file which
% contains all the details of GPIS model
% the second main output of this function is "cov_id" matrix which
% bascically tells which covariance function was optimum for training GPIS
% model

T_vector2 = input_cmg_mat;
[ererw,runs] = size(input_cmg_mat);

[lmv, gmv] = lmv_gmv_func(out_data,tar_data);

comb_inp_init{1,1} = lmv;  % LMV are the inputs of GP model

[comb_inp{1,1},inp_ind] = removeconstantrows(comb_inp_init{1,1});


u=1;

P1 = cell(u,1);
p1rows = cell(u,1);
P1n = cell(u,1);
p1s = cell(u,1);

%% normalize input and output data to [0 1]. Further, center the data to mean = 0;
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

well_num = 1;

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

dsfd = Target_cell;
Target_cell{1,1} = inpt_vect;
inpt_vect = dsfd{1,1};


    %% Divide data into training and testing set
    for data_div = 1 : 10
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
    
    [runs,num_dim] = size(y);   %y: parameters
    [runs,num_lvm] = size(x);  %x: lvms
    
    %% Train GPFS model for each output seperately
    if bUseGPLVM4GPIS == 0
         for jj = 1 : num_dim
            [hyp2{jj,1}, mu{jj,1}, s2{jj,1}, datastruc{jj,1}, cov_id{jj,1}] = train_gp_mult(x,y(:,jj),xs);
            sdev(jj,:) = 2.*sqrt(s2{jj,1});
            mu_store(jj,:) = mu{jj,1}';
         end 
    else
        %Using GPLVM based GPLVMIS to replace GPIS, which should be more
        %efficient especially when the number of uncertain parameters
        %(num_dim) is large. Actrually, here we can further discard the loop,
        %and only one gradients based optimization is required to get the
        %inverse sample, but to smoothly revise the code we still use the
        %loop, where num_lvm is typically much smaller than num_dim.
       %% Set up the GPLVM model(x: lvms; y: parameters)
        for jj = 1 : num_lvm
    %         [hyp2{jj,1}, mu{jj,1}, s2{jj,1}, datastruc{jj,1}, cov_id{jj,1}] = train_gp_mult(x,y(:,jj),xs);

            options = gpOptions('ftc');
            options.optimiser = 'optimiMinimize';%scg
            % options.kern = {'rbfard','white'};
            model = gpCreate(num_dim, 1, y, x(:,jj), options);
            model = gpOptimise(model, 1, -200);
            [mu{jj,1}, s2{jj,1}] = gpPosteriorMeanVar(model, ys);
            [hyp2{jj,1}, cov_id{jj,1}] = kernExtractParam(model.kern);
            datastruc{jj,1} = model.kern;        
            sdev(jj,:) = 2.*sqrt(s2{jj,1});
            mu_store(jj,:) = mu{jj,1}';
        end    
    end
    
    lb_test = mu_store - sdev;
    ub_test = mu_store + sdev;

    %% plot regression graph on the test set
% figure(11110)
% plotregression(T_test,mu_store,'GPIS Model Performance');
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



save inverse_model.mat



    
    
 
    
    
    
end