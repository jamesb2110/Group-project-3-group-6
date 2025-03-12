clear; clc;

a = 0.95;         
epsilon_range = linspace(0.01, 0.5, 100);
T_final = 1000;    
tspan = linspace(0, T_final, 10000);  
Y0 = [0; 0];    

x_vals = [];
c_vals = [];

for epsilon_val = epsilon_range
    model = @(t, Y) [
        (Y(1) - (Y(1)^3)/3 - Y(2)) / epsilon_val;  
        Y(1) + a;                              
    ];

    [t, Y] = ode45(model, tspan, Y0);
    
    X_vals = Y(:, 1);  
    
    last_points = X_vals(end-499:end);
    x_vals = [x_vals; last_points];
    c_vals = [c_vals; repmat(epsilon_val, 500, 1)];
end

figure;
plot(c_vals, x_vals, '.', 'MarkerSize', 5);
xlabel('epsilon');
ylabel('x');
grid on;
