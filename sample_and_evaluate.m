% muestreo perturbaciones, simulo K candidatos y calculo costo
function [S, E] = sample_and_evaluate(x0, lambda, U_nom, sigma, K, A, B, Cr, Yref, u_prev, c1, c2, u_max)
[nu, N] = size(U_nom);
qw = [3; 1];                              

% muestreo perturbaciones
E = sigma .* randn(nu, N, K); % ruido [w]

%claude additions, mejora ligeramente pero amuenta el overshoot. Mejora
%muy notable en ks mas bajos, pero al aumentar K el overshoot es mas
%grande a comparacion de las figuras sin el fix
half = floor(K/2);
E(:,:,half+1:2*half) = -E(:,:,1:half);  
E(:,:,1) = 0;                           

V = U_nom + E; % controles candidatos [w]
V = max(min(V,u_max),-u_max);                  
E = V - U_nom;                            

% simulacion y calculo de costo
Xk     = repmat(x0, 1, K);              
u_last = repmat(u_prev, 1, K);          
C1 = zeros(1, K); % costo de seguimiento
C2 = zeros(1, K); % costo de cambio de control

Sinv = 1 ./ sigma.^2;          
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