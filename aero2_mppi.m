% REFERNECIAS
% [mw]  mathworks, "local path planning using model predictive path integral"
% [w]   williams et al. 2017, "model predictive path integral control:
%       from theory to parallel computation"

clear; clc; close all

Dt  = 0.16743;
Mb  = 1.07;
Dm  = 2.4e-3;
Jp  = 23188500.45 / 1000 * 0.001^2;
Jy  = 23810415.95 / 1000 * 0.001^2;
g   = 9.81;
Ksp = 0.00744;
Dp  = 0.00199;
Kpp = 0.0025;
Kpy = 0.0015;
Dy  = 0.00192;
Kyy = 0.0068;
Kyp = -0.0033;
uno = -Ksp/Jp;   
dos = 0;  
tres = -Dp/Jp;   
cuatro = 0;
cinco = 0;       
seis = 0;  
siete = 0;       
ocho = -Dy/Jy;
nueve = (Kpp*Dt)/Jp;   
diez = (Kpy*Dt)/Jp;
once  = (Kyp*Dt)/Jy;  
doce = (Kyy*Dt)/Jy;
A_t = [0,0,1,0; 0,0,0,1; uno,dos,tres,cuatro; cinco,seis,siete,ocho];
B_t = [0 0; 0 0; nueve diez; once doce];
Cr  = [1 0 0 0; 0 1 0 0];
D_t = zeros(2,2);
Ts  = 0.1;
sys_d = c2d(ss(A_t, B_t, Cr, D_t), Ts, 'zoh');   
A_t = sys_d.A;   B_t = sys_d.B;   Cr = sys_d.C;
x0  = zeros(4,1);                               

c1 = 1; % peso del error de seguimiento
c2 = 1e-4;% peso del la suavidad cambio de control 

%% ==========
K_list = [10 1000];    
N_list = [10];    
sigma  = [3; 3];   % desviacion del ruido que le meto al control
lambda = 0.002;    % temperatura
T      = 210;      % pasos de simulacion
N_iter = 3;        % pasadas de mppi por paso (lo agregue pq no alcanzaba los peaks)
%% ======

for N = N_list
for K = K_list
    rng(0);                                         
    tt    = 0:T+N;
    Yfull = [0.3*sin(0.05*tt); 0.5*cos(0.05*tt)]; % referencia completa
    x = x0;  U = zeros(2,N);  u_prev = [0;0];
    Ylog = zeros(2,T);                               
    tic
    for t = 1:T
        Yref = Yfull(:, t:t+N); % ventana de referencia de este paso

        for it = 1:N_iter
            % muestreo K rollouts, los simulo y calculo el costo S
            [S,E] = sample_and_evaluate(x,lambda, U, sigma, K, A_t, B_t, Cr, Yref, u_prev, c1, c2);
            % saco los pesos con S y actualizo la secuencia U
            [U,w] = mppi_update(U, E, S, lambda);
        end

        % aplico el primer control
        u = U(:,1);
        x = step_dynamics(x, u, A_t, B_t);           
        U = [U(:,2:end), U(:,end)];                 
        u_prev = u;
        Ylog(:,t) = Cr*x;                            
    end

    t_ciclo = toc / T;                              
    rmse = sqrt(mean((Ylog - Yfull(:,2:T+1)).^2, 'all'));
    fprintf('N=%2d  K=%5d  t/ciclo=%7.2f ms  RMSE=%.4f\n', N, K, 1000*t_ciclo, rmse);

    figure;
    plot(1:T, Ylog', 1:T, Yfull(:,2:T+1)', '--');
    title(sprintf('N=%d, K=%d', N, K));
    xlabel('Time step'); ylabel('Output value');                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    legend('y_1','y_2','y_1 ref','y_2 ref');
    grid on;
end
end

% steps
function x_next = step_dynamics(x, u, A, B)
    x_next = A*x + B*u; %[mw]
end

% muestreo perturbaciones, simulo K candidatos y calculo costo
function [S, E] = sample_and_evaluate(x0, lambda, U_nom, sigma, K, A, B, Cr, Yref, u_prev, c1, c2)
    [nu, N] = size(U_nom);
    qw = [3; 1];                              

    % muestreo perturbaciones
    E = sigma .* randn(nu, N, K); % ruido [w]

    %claude additions, mejora ligeramente pero amuenta el overshoot. Mejora
    %muy notable en ks mas bajos, pero al aumentar K el overshoot es mas
    %grande a comparacion de las figuras sin el fix
    half = floor(K/2);
    E(:,:,half+1:2*half) = -E(:,:,1:half);   % antithetic pairs
    E(:,:,1) = 0;                            % keep the nominal as a candidate

    V = U_nom + E; % controles candidatos [w]
    V = max(min(V,24),-24);                  
    E = V - U_nom;                            

    % simulacion y calculo de costo
    Xk     = repmat(x0, 1, K);              
    u_last = repmat(u_prev, 1, K);          
    C1 = zeros(1, K); % costo de seguimiento
    C2 = zeros(1, K); % costo de cambio de control

    Sinv = 1 ./ sigma.^2;          % Sigma^-1 (diagonal)
    C3 = zeros(1, K);

    for t = 1:N
        Et = reshape(E(:,t,:), nu, K);
        C3 = C3 + lambda * sum((U_nom(:,t) .* Sinv) .* Et, 1); %It penalizes perturbations that push further in the direction the nominal control already points [w]
        Vt = reshape(V(:,t,:), nu, K);       
        Xk = A*Xk + B*Vt; % dinamica [w][mw]
        e  = Cr*Xk - Yref(:,t+1); % error de seguimiento [mw]
        du = Vt - u_last; % cambio de control [mw]
        C1 = C1 + sum(qw .* e.^2, 1); %[mw]
        C2 = C2 + sum(du.^2, 1); % [mw]
        u_last = Vt;
    end
    %S = (c1*C1 + c2*C2)'; % costo por rollout[mw]
    S = (c1*C1 + c2*C2 + C3)'; % costo por rollout[mw]

end

% Se calcula el peso en base al costo, se saca el promedio ponderado de todos los rollouts
function [U_new, w] = mppi_update(U_nom, E, S, lambda)
    K = numel(S);
    w = exp(-(S - min(S)) / lambda);% costo bajo = peso alto [w]
    w = w / sum(w);% normalizo [w]
    U_new = U_nom + sum(E .* reshape(w, 1, 1, K), 3); % promedio ponderado[w]
    U_new = max(min(U_new,24),-24); % saturacion
  
end