
clear; clc;


epsilon = 0.2;   
tau = 1.5;        
a = 0.5;         

c_range = linspace(-10.0, 10.0, 100);  
T_final = 300;    
t_transient = 200; 
t_span = [0, T_final]; 

history = @(t) [0.1, 0.1]; 

x_vals = cell(length(c_range), 1);
c_vals = cell(length(c_range), 1);

fprintf('开始计算耦合强度参数c的分岔图...\n');
progress_bar = waitbar(0, '计算中...');

for i = 1:length(c_range)
    c_val = c_range(i);
    
    waitbar(i/length(c_range), progress_bar, sprintf('计算 c = %.2f (%d/%d)', c_val, i, length(c_range)));

    model = @(t, Y, Z) [
        (Y(1) - (Y(1)^3) / 3 - Y(2) + c_val * (Z(1) - Y(1))) / epsilon;  
        Y(1) + a; 
    ];

    options = ddeset('RelTol', 1e-3, 'AbsTol', 1e-5);
    
    sol = dde23(model, tau, history, t_span, options);
    
    t_plot = linspace(t_transient, T_final, 150);
    Y_plot = deval(sol, t_plot);

    x_vals{i} = Y_plot(1, :)';
    c_vals{i} = repmat(c_val, length(t_plot), 1);
end

close(progress_bar);

x_all = vertcat(x_vals{:});
c_all = vertcat(c_vals{:});

figure;
plot(c_all, x_all, '.m', 'MarkerSize', 1);
xlabel('c');
ylabel('X');

grid on;
saveas(gcf, 'bifurcation_c.png');

