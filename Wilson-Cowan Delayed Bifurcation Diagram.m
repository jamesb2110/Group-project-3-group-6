clear;                           % clear variables
format compact
close all;                       % close figures

% Define the system right-hand side (RHS) and time constants
neuron_sys_rhs = @(xx, par) [ 
    -xx(1,1) + 1 / (1 + exp(-par(7) * (par(1) * xx(1,2) + par(2) * xx(2,3) + par(5)))); 
    -xx(2,1) + 1 / (1 + exp(-par(7) * (par(3) * xx(1,3) + par(4) * xx(2,2) + par(6))))];

ind_tau1=8;  % used later for continuation
ind_tau2=9; % used later for continuation

% Define the time constants (for DDEs)
neuron_tau = @() [8, 9];  % Time constants

% Set up the function handles for the DDE-Biftool
funcs = set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', neuron_tau);

% Initial guess for the steady state
stst.kind = 'stst';
stst.parameter = [-1, -0.4, -1, 0, 0.65, 0.5, 60, 1, 1];  % Start with par(7) = 60
stst.x = [0.5105467236749; 0.3468751427646];  % Adjust initial guess (you might need to tweak this)

% Set the continuation method (standard for DDEs)
method = df_mthod(funcs, 'stst');  % Use df_mthod for steady state continuation
method.stability.minimal_real_part = -1;

% Increase precision in the correction method
tol = 1e-8;  % Set tighter tolerance for p_correc
[stst, success] = p_correc(funcs, stst, [], [], method.point, tol);  % Correct the point for steady state

% Compute the stability of the steady state
stst.stability = p_stabil(funcs, stst, method.stability);

% Plot the stability
figure(1); clf;
p_splot(stst);  % Plot the stability diagram

method.stability.minimal_real_part=-2;
stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:
figure(2); clf;
p_splot(stst); % replot stability

% get an empty branch with ind_tau1 as a free parameter:
branch1=df_brnch(funcs,ind_tau1,'stst')
branch1.parameter
branch1.parameter.min_bound
% set bounds for continuation parameter
branch1.parameter.min_bound(1,:)=[ind_tau1 0];
branch1.parameter.max_bound(1,:)=[ind_tau1 1.5];
branch1.parameter.max_step(1,:)=[ind_tau1 0.02];
% use stst as a first branch point:
branch1.point=stst;

stst.parameter(ind_tau1)=stst.parameter(ind_tau1)+0.001;
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% use as a second branch point:
branch1.point(2)=stst;
branch1.method.continuation.plot=0;

% continue in one direction:
[branch1,s,f,r]=br_contn(funcs,branch1,100)
% turn the branch around:
branch1=br_rvers(branch1);
% continue in the other direction:
[branch1,s,f,r]=br_contn(funcs,branch1,100)

branch1.method.stability.minimal_real_part=-2;
branch1=br_stabl(funcs,branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,branch1)
figure(3); clf;
br_plot(branch1,xm,ym,'b'); % plot stability along branch:
ym.subfield='l0';
br_plot(branch1,xm,ym,'c');
plot([0 5],[0 0],'-.');
axis([0 5 -2 1.5]);
xlabel('tau1');ylabel('\Re\lambda');
% plot stability versus point number:
figure(4); clf;
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');

ind_hopf=find(arrayfun(@(x)real(x.stability.l0(1))>0,branch1.point),1,'last');
hopf=p_tohopf(funcs,branch1.point(ind_hopf));
method=df_mthod(funcs,'hopf'); % get hopf calculation method parameters:
method.stability.minimal_real_part=-1;
[hopf,success]=p_correc(funcs,hopf,ind_tau1,[],method.point) % correct hopf
first_hopf=hopf;                    % store hopf point in other variable for later use
hopf.stability=p_stabil(funcs,hopf,method.stability); % compute stability of hopf point
figure(5); clf;
p_splot(hopf);                     % plot stability of hopf point

branch2=df_brnch(funcs,[ind_tau1,ind_tau2],'hopf'); % use hopf point as first point of hopf branch:
branch2.parameter.min_bound(1,:)=[ind_tau1 0];
branch2.parameter.max_bound(1:2,:)=[[ind_tau1 1]' [ind_tau2 2]']';
branch2.parameter.max_step(1:2,:)=[[ind_tau1 0.02]' [ind_tau2 0.05]']';
branch2.point=hopf;
branch2.method.continuation.plot_progress=0;

hopf.parameter(ind_tau2)=hopf.parameter(ind_tau2)+0.001; % perturb hopf point
[hopf,success]=p_correc(funcs,hopf,ind_tau1,[],method.point); % correct hopf point, recompute stability
branch2.point(2)=hopf;                                 % use as second point of hopf branch:
branch2.method.continuation.plot_progress=0;
figure(6); clf;
[branch2,s,f,r]=br_contn(funcs,branch2,40);            % continue with plotting hopf branch:
branch2=br_rvers(branch2);                             % reverse Hopf branch
[branch2,s,f,r]=br_contn(funcs,branch2,30);            % continue in other direction
xlabel('tau1');ylabel('tau2');

branch3=df_brnch(funcs,[ind_tau1,ind_tau2],'hopf');
branch3.parameter=branch2.parameter;
branch3.point=hopf;
% perturb and correct:
hopf.parameter(ind_tau1)=hopf.parameter(ind_tau1)-0.05;
method.point.print_residual_info=0;
format short;
[hopf,success]=p_correc(funcs,hopf,ind_tau2,[],method.point);
branch3.point(2)=hopf; % use as second branch point:
% continue branch of hopf points on two sides:
branch3.method.continuation.plot_progress=0;
figure(6); clf;

% Perform continuation in both directions and capture branch points
[branch3, s, f, r] = br_contn(funcs, branch3, 1000);

% Reverse branch to go in the opposite direction
branch3 = br_rvers(branch3);
[branch3, s, f, r] = br_contn(funcs, branch3, 1000);

% Extract tau1 and tau2 values from the branch points
tau1_values = [];
tau2_values = [];

% Loop through each point in the branch and extract the parameter values for tau1 and tau2
for i = 1:length(branch3.point)
    tau1_values = [tau1_values; branch3.point(i).parameter(ind_tau1)];
    tau2_values = [tau2_values; branch3.point(i).parameter(ind_tau2)];
end

% Check if tau1_values and tau2_values are the same length
if length(tau1_values) == length(tau2_values)
    % Plot the clean bifurcation line for Hopf points
    plot(tau1_values, tau2_values, 'black', 'LineWidth', 2);
    
    % Set axis limits to focus on the range of interest
    axis([0 1 0 1]);  % Adjust as needed for the range of tau1 and tau2

    % Add labels and title for clarity
    ax=gca;
    ax.FontSize = 20
    xlabel('$\tau_1$', 'Interpreter', 'latex','FontSize',25);
    ylabel('$\tau_2$', 'Interpreter', 'latex','FontSize',25);

else
    error('Tau1 and Tau2 arrays have different lengths. Check the branch points structure.');
end


