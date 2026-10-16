% Se calcula el peso en base al costo, se saca el promedio ponderado de todos los rollouts
function [U_new, w] = mppi_update(U_nom, E, S, lambda, u_max)
K = numel(S);
w = exp(-(S - min(S)) / lambda);% costo bajo = peso alto [w]
w = w / sum(w);% normalizo [w]
U_new = U_nom + sum(E .* reshape(w, 1, 1, K), 3); % promedio ponderado[w]
U_new = max(min(U_new,u_max),-u_max); 
end