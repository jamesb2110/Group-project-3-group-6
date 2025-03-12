
clear; clc;

epsilon = 0.1;   
tau = 1.0;       
c = -0.5;         

a_range = linspace(-2.0, 2.0, 50);  
T_final = 300;    
t_transient = 200;
t_span = [0, T_final]; 

history = @(t) [0.1, 0.1];  

x_vals = cell(length(a_range), 1);
a_vals = cell(length(a_range), 1);

fprintf('开始计算参数a的分岔图...\n');
progress_bar = waitbar(0, '计算中...');

for i = 1:length(a_range)
    a_val = a_range(i);

    waitbar(i/length(a_range), progress_bar, sprintf('计算 a = %.2f (%d/%d)', a_val, i, length(a_range)));

    model = @(t, Y, Z) [
        (Y(1) - (Y(1)^3) / 3 - Y(2) + c * (Z(1) - Y(1))) / epsilon; 
        Y(1) + a_val; 
    ];

    options = ddeset('RelTol', 1e-3, 'AbsTol', 1e-5);

    sol = dde23(model, tau, history, t_span, options);

    t_plot = linspace(t_transient, T_final, 150);
    Y_plot = deval(sol, t_plot);

    x_vals{i} = Y_plot(1, :)';
    a_vals{i} = repmat(a_val, length(t_plot), 1);
end

close(progress_bar);

x_all = vertcat(x_vals{:});
a_all = vertcat(a_vals{:});

figure;
plot(a_all, x_all, '.k', 'MarkerSize', 1);
xlabel('a');
ylabel('X');
grid on;

saveas(gcf, 'bifurcation_a.png');

