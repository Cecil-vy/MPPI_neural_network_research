function x_next = step_dynamics(x, u, A, B)
x_next = A*x + B*u; %[mw]
end