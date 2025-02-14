%% Results

clear
clc

load Hearthquakes_Forecasts

MSE_p = (sum(1/size(y_t_all(WE+1:end,:),1)*(y_t_all(WE+1:end,:)...
    -lambda_t_p_forecast(WE:end,:)).^2))
MSE_par = (sum(1/size(y_t_all(WE+1:end,:),1)*(y_t_all(WE+1:end,:)...
    -lambda_t_par_forecast(WE:end,:)).^2))
MSE_partrend = (sum(1/size(y_t_all(WE+1:end,:),1)*(y_t_all(WE+1:end,:)...
    -lambda_t_partrend_forecast(WE:end,:)).^2))
MSE_tpar = (sum(1/size(y_t_all(WE+1:end,:),1)*(y_t_all(WE+1:end,:)...
    -lambda_t_tpar_forecast(WE:end,:)).^2))

FS_p = sum(1/size(y_t_all(WE+1:end,1),1)*(y_t_all(WE+1:end,1).*log(lambda_t_p_forecast(WE:end,:))...
    -(lambda_t_p_forecast(WE:end,:))).^2)
FS_par = sum(1/size(y_t_all(WE+1:end,1),1)*(y_t_all(WE+1:end,1).*log(lambda_t_par_forecast(WE:end,:))...
    -(lambda_t_par_forecast(WE:end,:))).^2)
FS_partrend = sum(1/size(y_t_all(WE+1:end,1),1)*(y_t_all(WE+1:end,1).*log(lambda_t_partrend_forecast(WE:end,:))...
    -(lambda_t_partrend_forecast(WE:end,:))).^2)
FS_tpar = sum(1/size(y_t_all(WE+1:end,1),1)*(y_t_all(WE+1:end,1).*log(lambda_t_tpar_forecast(WE:end,:))...
    -(lambda_t_tpar_forecast(WE:end,:))).^2)

TickXAxis = [12:24:size(lambda_t_p(1200:end,1),1)];

figure(3)
% subplot(2,1,1)
p1 =  plot(y_t_all(1200:end,1), 'LineWidth',LineWidth,'Color', 'b' ,'LineStyle', '-')
hold on
p2 = plot(lambda_t_p_forecast(1200:end,1), 'LineWidth',LineWidth,'LineStyle', '--')
p3 = plot(lambda_t_par_forecast(1200:end,1), 'LineWidth',LineWidth,'LineStyle', '--')
p4 = plot(lambda_t_partrend_forecast(1200:end,1), 'LineWidth',LineWidth ,'LineStyle', '--')
p5 = plot(lambda_t_tpar_forecast(1200:end,1), 'LineWidth',LineWidth,'LineStyle', '--')
title('\textbf{Monthly number of hearthquakes in Italy with forecasts}','interpreter', 'latex')
% leg1 = legend('Tornadoes ($y_t$)');
leg1 = legend('$y_t$','$\widehat{\lambda}$ - Poisson','$\widehat{\lambda}_{T+1}$ - PAR','$\widehat{\lambda}_{T+1}$ - PARX','$\widehat{\lambda}_{T+1}$ - SETPAR');
set(gca,'FontSize',FontSizeGraph);
set(leg1,'Interpreter','latex');
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',{'1/2001','1/2003','1/2005','1/2007','1/2009','1/2011', ...
                      '1/2013','1/2015','1/2017','1/2019','1/2021','1/2023'});  
set(gca,'XTick',TickXAxis );
leg1.FontSize = FontSize_Legend;   
axis tight


