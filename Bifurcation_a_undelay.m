clear; clc;

epsilon = 1;  
c = 1;           
a_range = linspace(0.5, 1.5, 100);
T_final = 100000;     
tspan = linspace(0, T_final, 1000000);  
Y0 = [0; 0];    

x_vals = [];
c_vals = [];

for a_val = a_range

    model = @(t, Y) [
        (Y(1) - (Y(1)^3)/3 - Y(2)) / epsilon;  
        Y(1) + a_val;                              
    ];

    [t, Y] = ode45(model, tspan, Y0);
    
    X_vals = Y(:, 1); 
    
    last_points = X_vals(end-499:end);
    x_vals = [x_vals; last_points];
    c_vals = [c_vals; repmat(a_val, 500, 1)];
end

figure;
plot(c_vals, x_vals, '.', 'MarkerSize', 5);
xlabel('a');
ylabel('x');
grid on;
