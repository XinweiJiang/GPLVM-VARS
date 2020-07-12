function [hyp2, mu, s2, datastruc] = train_gp_mult_constcov(x,y,xs,covid)

% this function is same as train_gp_mult.m. However, this function has one
% additinal parameter covid, which means the index of a particular
% covariane function that was mostly optimum for last 10 iterations of
% GP-VARS.


[run,D] = size(x);
[run,num_dim] = size(y);
%      mean = {@meanSum,{@meanLinear,@meanConst}}; a = 1/5; b = 1;               % empty: don't use a mean function
%   covfunc ={@covSEiso};              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood

L = rand(D,1); sf = rand();
ell = rand(); sf = rand();
al = 2;

cgi = {'covSEiso'};  hypgi = log([ell;sf]);    % isotropic Gaussian
ccard = {'covMaternard',3}; L = rand(D,1); sf = rand(); hypard = log([L;sf]);
cciso = {'covMaterniso',3}; ell = rand(); sf = rand();  hypiso = log([ell;sf]);

% cov = {'covMaternard',3};                             % setup the GP


for mm = covid : covid
    L = rand(D,1);
    switch mm        case 1
            cov = {'covSEisoU'};  hypcov = log([ell]);    % isotropic Gaussian
        case 2
            cov = {'covMaternard',3};  hypcov = log([L;sf]);
        case 3
            cov = {'covMaterniso',3};  hypcov = log([ell;sf]);
        case 4
            cov = {@covNNone}; hypcov = log([ell;sf]);           % neural network
        case 5
            cov = {'covRQiso'};  hypcov = log([ell;sf;al]); % ration. quad.
        case 6
            cov1 = {'covMaterniso',3}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 7
            cov1 = {'covMaternard',3}; hypcov1 = log([L;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 8
            cov1 = {@covNNone}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covProd',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
        case 9
            cov1 = {@covNNone}; hypcov1 = log([ell;sf]);
            cov2 = {'covSEiso'}; hypcov2 = log([ell;sf]);
            cov = {'covSum',{cov1,cov2}};  hypcov = [hypcov1;hypcov2];
            
    end
    
    covcell{mm,1} = cov;
    hypcov_store = hypcov;
    
    for iiv = 1 : 3
        %     L = rand(D,1).*iiv ;sf = rand()*iiv ; ell = abs(rand()*iiv - 0.5);                             % setup the GP
        %     hyp0.cov  = log([L;sf]);
        hyp0.cov = hypcov_store.* abs((rand()*iiv) - 0.5);
        % hyp0.cov  = [];
        % hyp0.cov  = [0;0];
        %    mean = {@meanSum,{@meanLinear,@meanConst}}; a = rand(D,1); b = 0.5;       % m(x) = a*x+b
        mean = [];
        %   hyp0.mean = [a;b];
        hyp0.mean = [];
        
        hyp0.lik  = -1;
        covfunc = cov;
        
        %   hyp = struct('cov', [0;0],'mean', [a;b], 'lik', -1);
        [hyp2_cell{iiv,mm},~,fval(iiv,mm)] = minimize(hyp0, @gp, -200, @infGaussLik , mean, covfunc, likfunc, x, y);
        
    end
end

[min_fval_row,min_fval_loc] = min(fval);
[min_fval_col,min_fval_loc_col] = min(min_fval_row);

col_num = covid;

row_num = min_fval_loc(1,col_num);


hyp2 = hyp2_cell{row_num,col_num};

covfunc = covcell{col_num,1};

datastruc.inf = {'infGaussLik'};
datastruc.mean = mean;
datastruc.covfunc = covfunc;
datastruc.likfunc = likfunc;

ep = 0.001;
err_msg = 'error';
stat_chng = 0;
[ui,vi] = size(xs);
while strcmp(err_msg,'error') == 1
    try
        [mu s2] = gp(hyp2, @infGaussLik , mean, covfunc, likfunc, x, y, xs);
        err_msg = 'good';
    catch nosense
        err_msg = 'error';
        xs_t = xs';
        
       xs_t = xs_t + ones(vi,ui).*ep;
       xs = xs_t';
        stat_chng = stat_chng+1;
    end
end

end