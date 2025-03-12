clear; clc;

epsilon = 0.1;   
a = 0.5;          
c = -1.0;         


tau_range = linspace(0.01, 4.0, 50);  
T_final = 300;    
t_transient = 200; 
t_span = [0, T_final]; 


history = @(t) [0.1, 0.1];  


x_vals = cell(length(tau_range), 1);
tau_vals = cell(length(tau_range), 1);


fprintf('开始计算延迟参数τ的分岔图...\n');
progress_bar = waitbar(0, '计算中...');

for i = 1:length(tau_range)
    tau_val = tau_range(i);

    waitbar(i/length(tau_range), progress_bar, sprintf('计算 τ = %.2f (%d/%d)', tau_val, i, length(tau_range)));

    model = @(t, Y, Z) [
        (Y(1) - (Y(1)^3) / 3 - Y(2) + c * (Z(1) - Y(1))) / epsilon;  
        Y(1) + a;  
    ];


    options = ddeset('RelTol', 1e-3, 'AbsTol', 1e-5);
    

    sol = dde23(model, tau_val, history, t_span, options);
    
    t_plot = linspace(t_transient, T_final, 150);
    Y_plot = deval(sol, t_plot);
    

    x_vals{i} = Y_plot(1, :)';
    tau_vals{i} = repmat(tau_val, length(t_plot), 1);
end


close(progress_bar);

x_all = vertcat(x_vals{:});
tau_all = vertcat(tau_vals{:});

figure;
plot(tau_all, x_all, '.b', 'MarkerSize', 1);
xlabel('τ');
ylabel('X');
grid on;

saveas(gcf, 'bifurcation_tau.png');

