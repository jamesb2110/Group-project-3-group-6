clear; clc;

a = 0.7;          
tau = 1.0;       
c = -0.5;        

epsilon_log = linspace(log10(0.01), log10(1.0), 50);
epsilon_range = 10.^epsilon_log;  
T_final = 300;    
t_transient = 200; 
t_span = [0, T_final]; 


history = @(t) [0.1, 0.1];  


x_vals = cell(length(epsilon_range), 1);
epsilon_vals = cell(length(epsilon_range), 1);


fprintf('开始计算时间尺度参数ε的分岔图...\n');
progress_bar = waitbar(0, '计算中...');

for i = 1:length(epsilon_range)
    epsilon_val = epsilon_range(i);
    

    waitbar(i/length(epsilon_range), progress_bar, sprintf('计算 ε = %.4f (%d/%d)', epsilon_val, i, length(epsilon_range)));
    

    model = @(t, Y, Z) [
        (Y(1) - (Y(1)^3) / 3 - Y(2) + c * (Z(1) - Y(1))) / epsilon_val;  
        Y(1) + a;  
    ];

    
    rel_tol = min(1e-3, epsilon_val/10);
    abs_tol = min(1e-5, epsilon_val/100);
    options = ddeset('RelTol', rel_tol, 'AbsTol', abs_tol);
    

    sol = dde23(model, tau, history, t_span, options);
    

    t_plot = linspace(t_transient, T_final, 150);
    Y_plot = deval(sol, t_plot);

    x_vals{i} = Y_plot(1, :)';
    epsilon_vals{i} = repmat(epsilon_val, length(t_plot), 1);
end


close(progress_bar);


x_all = vertcat(x_vals{:});
epsilon_all = vertcat(epsilon_vals{:});


figure;
semilogx(epsilon_all, x_all, '.r', 'MarkerSize', 1);
xlabel('ε');
ylabel('X');

grid on;

saveas(gcf, 'bifurcation_epsilon.png');

