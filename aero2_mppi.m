% REFERNECIAS
% [mw]  mathworks, "local path planning using model predictive path integral"
% [w]   williams et al. 2017, "model predictive path integral control:
%       from theory to parallel computation"

clear; clc; 

modelo_aero2;   

u_max = 18;     

c1 = 1; 
c2 = 1e-4;

%% ==========
barrido = true;  
if barrido
    K_list = [2 5 10 25 50 100 250 1000];
    N_list = [5 10 20];
else
    K_list = 1000;
    N_list = 20;
end

seeds  = 1:10;
sigma  = [3; 3];   
lambda = 0.002;    
T   = 100;  
N_iter = 3;        % pasadas de mppi por paso (lo agregue pq no alcanzaba los peaks)


if ~exist('resultados','dir'), mkdir('resultados'); end

res = [];   

for N = N_list
    for K = K_list
        tt    = 0:T+N;
        tsec  = tt*Ts;
        Yfull = [interp1([0 2 4 6 8], [3 -4.5 6 -0.5 4], tsec, 'previous', 'extrap'); ...
                 interp1([0 3 6 9],   [-5 8 -2 5],       tsec, 'previous', 'extrap')] * pi/180;
        rmse    = zeros(1, numel(seeds));
        tc      = zeros(1, numel(seeds));
        rmse_ch = zeros(2, numel(seeds));  
        upk     = zeros(1, numel(seeds));   
        sat     = zeros(1, numel(seeds));  

        for s = 1:numel(seeds)
            rng(seeds(s));
            x = x0;  U = zeros(2,N);  u_prev = [0;0];
            Ylog = zeros(2,T);
            Ulog = zeros(2,T);
            tic
            for t = 1:T
                Yref = Yfull(:, t:t+N);

                for it = 1:N_iter
                    [S,E] = sample_and_evaluate(x,lambda, U, sigma, K, A_t, B_t, Cr, Yref, u_prev, c1, c2, u_max);
                    [U,w] = mppi_update(U, E, S, lambda, u_max);
                end

                u = U(:,1);
                x = step_dynamics(x, u, A_t, B_t);
                U = [U(:,2:end), U(:,end)];
                u_prev = u;
                Ylog(:,t) = Cr*x;
                Ulog(:,t) = u;
            end

            tc(s)   = toc / T;
            rmse(s) = sqrt(mean((Ylog - Yfull(:,2:T+1)).^2, 'all'));
            rmse_ch(:,s) = sqrt(mean((Ylog - Yfull(:,2:T+1)).^2, 2)) * 180/pi;
            upk(s) = max(abs(Ulog), [], 'all');
            sat(s) = mean(abs(Ulog) >= u_max, 'all') * 100;
        end

        res(end+1,:) = [N, K, mean(rmse), std(rmse), 1000*mean(tc), mean(rmse_ch,2)', mean(upk), mean(sat)];
        fprintf('N=%2d  K=%5d  t/ciclo=%7.2f ms  RMSE cab=%.2f deg  gui=%.2f deg  V pico=%.1f  sat=%.1f%%\n', N, K, 1000*mean(tc), mean(rmse_ch(1,:)), mean(rmse_ch(2,:)), mean(upk), mean(sat));

        if N==20 && K==1000
            ts = (1:T)*Ts;
            figure;
            subplot(3,1,1); plot(ts, Yfull(1,2:T+1)*180/pi, 'k--', ts, Ylog(1,:)*180/pi, 'b');
            ylabel('Cabeceo [deg]'); legend('referencia','MPPI'); grid on
            title(sprintf('Cabeceo | error RMS: MPPI %.2f deg', rmse_ch(1,end)))
            subplot(3,1,2); plot(ts, Yfull(2,2:T+1)*180/pi, 'k--', ts, Ylog(2,:)*180/pi, 'b');
            ylabel('Guiñada [deg]'); grid on
            title(sprintf('Guiñada | error RMS: MPPI %.2f deg', rmse_ch(2,end)))
            subplot(3,1,3); plot(ts, Ulog'); yline(u_max,':'); yline(-u_max,':');
            ylabel('Comando [V]'); xlabel('tiempo [s]'); legend('rotor de cabeceo','rotor de guiñada'); grid on
            title(sprintf('Comando | pico %.1f V de %d V disponibles', upk(end), u_max))

            fprintf('error RMS cabeceo : MPPI %.2f deg\n', rmse_ch(1,end));
            fprintf('error RMS guinada : MPPI %.2f deg\n', rmse_ch(2,end));
            fprintf('comando pico : MPPI %.1f V\n', upk(end));
            fprintf('saturacion : %.1f %% de las muestras\n', sat(end));
        end
    end
end

tabla = array2table(res, 'VariableNames', {'N','K','RMSE_mean','RMSE_std','ms_ciclo','RMSE_cab_deg','RMSE_gui_deg','V_pico','sat_pct'})
writetable(tabla, fullfile('resultados','baseline_mppi.csv'));