function [lambda_t] = Filter_TPAR_FORECAST(theta, y_tr, r)

%%% Get Tt
Tt = size(y_tr,1);

%%% Set Filters
log_lambda_t = zeros(Tt,1);
lambda_t = zeros(Tt,1);

%%% Parameters
omega_1 = theta(1);
beta_1 = theta(2);
alpha_1 = theta(3);
omega_2 = theta(4);
beta_2 = theta(5);
alpha_2 = theta(6);

%%% Initial values
% log_lambda_t(1) = omega/(1-beta);
lambda_t(1) = y_tr(1);

% if y(1) >= r
% log_lambda_t(1) = omega_1/(1-beta_1);
% lambda_t(1) = exp(log_lambda_t(1));
% else
% log_lambda_t(1) = omega_2/(1-beta_2);
% lambda_t(1) = exp(log_lambda_t(1));
% end

if y_tr(1) >= r
log_lambda_t(2) = omega_1/(1-beta_1);
lambda_t(2) = exp(log_lambda_t(2));
else
log_lambda_t(2) = omega_2/(1-beta_2);
lambda_t(2) = exp(log_lambda_t(2));
end

    %%% Updating
    for t = 2:(Tt)
        
        inno_t = (y_tr(t) - lambda_t(t))/(lambda_t(t)); 
        
        if y_tr(t-1)>=r
        log_lambda_t(t+1) = omega_1 + beta_1*log_lambda_t(t) + alpha_1*inno_t;
        else
        log_lambda_t(t+1) = omega_2 + beta_2*log_lambda_t(t) + alpha_2*inno_t;
        end
        
    lambda_t(t+1) = exp(log_lambda_t(t+1));
    end
end