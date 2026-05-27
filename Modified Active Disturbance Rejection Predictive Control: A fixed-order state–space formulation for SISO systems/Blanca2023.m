% Paper Simulation
% Modified Active Disturbance Rejection Predictive Control: A fixed-order state–space formulation for SISO systems
% https://doi.org/10.1016/j.isatra.2023.08.011

clear; clc; close all;
%% Uncertainties setup
uncertain = cell(3,2);
uncertain{1,1} = [0 0]; uncertain{1,2} = [-0.2 0.2];    % Simulation 1
uncertain{2,1} = [0]; uncertain{2,2} = [-0.2];          % Simulation 2
uncertain{3,1} = [0 0]; uncertain{3,2} = [0.2 -0.02];   % Simulation 3

Sim_case = cell(2);
Sim_case{1} = 'Nominal'; Sim_case{2} = 'With uncertainty';

%% Simulation setup
for simulation = 1:3
for aux = 1:2
switch simulation
    case 1  % Example 4.1
        % Model
        Gm_c = tf(2.5*(1+uncertain{1,aux}(1)),[0.9*(1+uncertain{1,aux}(2)) 1 0]);
        Ts = 0.05;
        % Constraints
        u_max = 24;
        u_min = -24;
        du_max = 5;
        du_min = -5;
        % ESO parameters
        tau = 0.9;
        b0 = 2.5/tau;
        omega_o = 20;
        % MPC parameters
        N = 40;
        Nu = 9;
        gamma = 1;
        lambda = 0.1;
        eps_1 = 1e5;
        eps_2 = 1e5;
        % Simulation parameters
        Tsim = 20;
        dist = 0;
        T_dist = 0;
        Fig_title = 'Example 4.1: DC motor';
    case 2  % Example 4.2.1
        % Model
        alpha = 0.5 * (1 + uncertain{2,aux});
        p1 = [1, 1]; p2 = [alpha, 1]; p3 = [alpha^2, 1]; p4 = [alpha^3, 1];
        den_c = conv(conv(p1, p2), conv(p3, p4));
        num_c = 1;
        Gm_c = tf(num_c, den_c);
        Ts = 0.1;
        % Constraints
        u_max = 1.2;
        u_min = -1.2;
        du_max = 0.5;
        du_min = -0.5;
        % ESO parameters
        tau = 1;
        b0 = 2.4;
        omega_o = 5;
        % MPC parameters
        N = 50;
        Nu = 5;
        gamma = 1;
        lambda = 0.1;
        eps_1 = 1e5;
        eps_2 = 1e5;
        % Simulation parameters
        Tsim = 20;
        dist = 1;
        T_dist = 10;
        Fig_title = 'Example 4.2.1: Fourth-Order Benchmark System';
    case 3  % Example 4.2.2
        % Model
        beta_val = 1 * (1 + uncertain{3,aux}(1)); 
        tau_p = 1 * (1 + uncertain{3,aux}(2));
        num_c = [-beta_val, 1];
        den_c = [tau_p^3, 3*tau_p^2, 3*tau_p, 1];
        Gm_c = tf(num_c, den_c);
        Ts = 0.1;
        % Constraints
        u_max = 5;
        u_min = -5;
        du_max = 1;
        du_min = -1;
        y_max = 2; 
        % ESO parameters
        tau = 1.62;
        b0 = 4.5;
        omega_o = 12;
        % MPC parameters
        N = 85;
        Nu = 20;
        gamma = 0.001;
        lambda = 5.5;
        eps_1 = 1e5;
        eps_2 = 1e5;
        % Simulation parameters
        Tsim = 30;
        dist = 1;
        T_dist = 15;
        Fig_title = 'Example 4.2.2: Third-Order Benchmark System';
end

% Discretization
Gm_d = c2d(Gm_c,Ts,'zoh');
[A_p,B_p,C_p,D_p] = ssdata(ss(Gm_d));

%% Disturbance rejector - Extended State Observer
z_o = exp(-omega_o*Ts);
a = exp(-Ts/tau);

A_o = [1, tau*(1 - a), tau*Ts - tau^2*(1 - a);
       0, a,           tau*(1 - a);
       0, 0,           1];

B_o = b0 * [tau*Ts - tau^2*(1 - a);
            tau*(1 - a);
            0];

C_o = [1, 0, 0];

lo1 = 1 - (z_o^3) / a;
lo3 = ((1 - z_o)^3) / (tau * Ts * (1 - a));
lo2 = (2*a - lo1*(1 + a) + lo3*(tau^2*(1 - a) - a*tau*Ts) - 3*z_o^2 + 1) / (tau*(1 - a));

L_o = [lo1; lo2; lo3];

%% MPC
A_m = A_o(1:2, 1:2);
B_m = B_o(1:2);
C_m = [1, 0];

P = zeros(N,2);
V_B = zeros(N,1);
G = zeros(N,Nu);

% Prediction Matrices
for i = 1:N
    P(i, :) = C_m * (A_m^i);
    
    sum_AB = zeros(2, 1);
    for j = 0:(i-1)
        sum_AB = sum_AB + (A_m^j) * B_m;
    end
    V_B(i, 1) = C_m * sum_AB;
    
    for j = 1:Nu
        if i >= j
            sum_G = zeros(2, 1);
            for k = 0:(i-j)
                sum_G = sum_G + (A_m^k) * B_m;
            end
            G(i, j) = C_m * sum_G;
        end
    end
end

H_du = 2 * (G' * gamma* eye(N) * G + lambda * eye(Nu));

% Adding the slack variables to cost function
n_eq = 2; % For a 2nd order system
P_eq = zeros(n_eq, 2);
V_B_eq = zeros(n_eq, 1);
G_eq = zeros(n_eq, Nu);

for idx = 1:n_eq
    i = N + idx;
    P_eq(idx, :) = C_m * (A_m^i);
    
    sum_AB = zeros(2, 1);
    for j = 0:(i-1)
        sum_AB = sum_AB + (A_m^j) * B_m;
    end
    V_B_eq(idx, 1) = C_m * sum_AB;
    
    for j = 1:Nu
        if i >= j
            sum_G = zeros(2, 1);
            for k = 0:(i-j)
                sum_G = sum_G + (A_m^k) * B_m;
            end
            G_eq(idx, j) = C_m * sum_G;
        end
    end
end
H_aug = blkdiag(H_du, 2 * eps_2 * eye(n_eq));

%% Simulation
display("Simulating the results from " + Fig_title + " - " + Sim_case{aux} + " case")

out = sim('simu_Blanca2023.slx');

% Figures
switch aux
    case 1
        figure(2*simulation-1)
        plot(out.t, out.ref, 'Color', '#A3A3A3', 'LineWidth', 1); hold on;
        plot(out.t, out.y, ':', 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        grid on;

        figure(2*simulation)
        subplot(2,1,1);
        plot(out.t, out.u, ':', 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        grid on;

        subplot(2,1,2);
        plot(out.t, out.du_o, ':', 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        grid on;
    case 2
        figure(2*simulation-1)
        plot(out.t, out.y, 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        title(Fig_title);
        xlabel('Time (s)'); ylabel('Output');
        legend('Reference','Nominal','w/ Uncertainty')
        grid on;

        figure(2*simulation)
        subplot(2,1,1);
        plot(out.t, out.u, 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        yline(u_max, '--k');
        yline(u_min, '--k');
        xlabel('Time (s)'); ylabel('Manipulated Variable');
        grid on;
        subplot(2,1,2);
        plot(out.t, out.du_o, 'Color', '#008B8B', 'LineWidth', 1.5); hold on;
        yline(du_max, '--k');
        yline(du_min, '--k');
        xlabel('Time (s)');
        ylabel('Rate of input');
        grid on;
        sgtitle(Fig_title)
end

end

end

display("Code author: Jose Sergio Cruz Dantas Junior - UFC")