%% This code replicates Section 4.2 of Angelini and Costantini (2025, ANOR)

clear all
close all
clc
format long

options = optimset('TolX', 1e-7, 'Maxiter', 50000, 'MaxFunEvals', 50000, 'HessUpdate', 'bfgs');

load Temperature_cities.txt

y_t_all = temperature_cities;

for j = 1 : 6
dummies(j+4:12:size(temperature_cities,1),j) = 1;
end

lambda_t_p_forecast = nan(size(y_t_all,1)-1,6);
lambda_t_par_forecast = nan(size(y_t_all,1)-1,6);
lambda_t_partrend_forecast = nan(size(y_t_all,1)-1,6);
lambda_t_tpar_forecast = nan(size(y_t_all,1)-1,6);

WE = 517;

for s = 1 : 6
    
for j = WE:size(y_t_all,1)-1

y_t = y_t_all(j-WE+1:j,s);
    
%% POISSON
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

lb_bench = [0];
ub_bench = [10];

theta_p = [2];

[theta_mle_p, loglik_mle_p, ~, ~, ~, ~, hessian_p] =...
            fmincon(@LogLik_P,theta_p,A,b,Aeq,beq,lb_bench,ub_bench,nonlcon,options,y_t);
               
log_lambda_t = theta_mle_p;
lambda_t = exp(log_lambda_t);   
lambda_t_p_forecast(j,s) = lambda_t;

%% PAR
lb_bench = [-1 0 -1];
ub_bench = [1 0.999 1];

theta_par = [0 0.8 0.5];

LogLik_PAR(theta_par,y_t)

[theta_mle_par, loglik_mle_par, ~, ~, ~, ~, hessian_par] =...
            fmincon(@LogLik_PAR,theta_par,A,b,Aeq,beq,lb_bench,ub_bench,nonlcon,options,y_t);
       
[lambda_t_par] = Filter_PAR_FORECAST(theta_mle_par, y_t);
lambda_t_par_forecast(j,s) = lambda_t_par(end);

%% PAR + trend
lb_bench = [-1 0 -1 0 ones(1,6)*-0.5];
ub_bench = [1 0.999 1 0.1 ones(1,6)*5];

theta_partrend = [theta_mle_par 0.000001 ones(1,6)*1];

linear_trend = (1:size(y_t,1)+1)';

LogLik_PARX(theta_partrend,y_t,[linear_trend, dummies(j-WE+1:j+1,:)])

[theta_mle_partrend, loglik_mle_partrend, ~, ~, ~, ~, hessian_partrend] =...
            fmincon(@LogLik_PARX,theta_partrend,A,b,Aeq,beq,lb_bench,ub_bench,nonlcon,options,y_t,[linear_trend, dummies(j-WE+1:j+1,:)]);
       
[lambda_t_partrend] = Filter_PARX_FORECAST(theta_mle_partrend, y_t,[linear_trend, dummies(j-WE+1:j+1,:)]);
lambda_t_partrend_forecast(j,s) = lambda_t_partrend(end);

%% TPAR
lb_bench = [-1 0 -1 -1 0 -1];
ub_bench = [1 0.999 1 1 0.999 1];

theta_tpar = [0 0.8 0.1 0 0.8 0.1];

for r = 1:5
[~, loglik_mle_tpar_iter(r,:), ~, ~, ~, ~, ~,] =...
            fmincon(@LogLik_TPAR,theta_tpar,A,b,Aeq,beq,lb_bench,ub_bench,nonlcon,options,y_t,r);
end

best_iter_r = find(loglik_mle_tpar_iter==min(loglik_mle_tpar_iter));

[theta_mle_tpar, loglik_mle_tpar(r,:), ~, ~, ~, ~, hessian_tpar] =...
            fmincon(@LogLik_TPAR,theta_tpar,A,b,Aeq,beq,lb_bench,ub_bench,nonlcon,options,y_t,best_iter_r);

[lambda_t_tpar] = Filter_TPAR_FORECAST(theta_mle_tpar, y_t,best_iter_r);
lambda_t_tpar_forecast(j,s) = lambda_t_tpar(end);

end

end

%% SAVE AND PRINT RESULTS
save('Temperature_Forecasts')

Results_Temperature_Forecasting