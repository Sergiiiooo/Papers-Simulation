% Paper simulation
% https://doi.org/10.1016/j.ifacol.2024.08.002

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% 2nd order model
% Continuous time
s = tf('s');
Ps = (-0.8*s+1)/(1.5*s+1)^2

% Sampling time
Ts = 0.1;

% Discrete time
Pz = c2d(Ps,Ts,'zoh')
[num,den] = tfdata(Pz,'v');
num = round(num, 3, 'decimals');
den = round(den, 3, 'decimals');
% Pz = tf(num,den,Ts)

%% Controller tuning parameters
% Prediction horizon
N = 20;
% Control horizon
Nu = 1;

% Weight values
lambda = 0;
lambda_e = 1000;
Q = eye(N);
for aux=1:8
    Q(aux,aux) = 1;
end

% Constraints
du_min = -0.5;
du_max = 0.5;
u_min = 0;
u_max = 0.9;
y_min = 0;
y_max = 0.7;

%% CARIMA model
Delta = tf([1 -1],[1],Ts,'Variable','z^-1'); % Δ = 1 - z^-1
[delta, ~] = tfdata(Delta, 'v');

% Polynomials b(z^-1) and a(z^-1)
[b,a] = tfdata(Pz,'v');

% Disturbance rejection polynomial c(z^-1) 
alf = 0.0;
c  = conv([1 -alf],[1 -alf]);

a_til = conv(delta,a); % a_til = Δ*a

Max = max([length(a_til)-1,length(b)-1,length(c)-1]);
a_til = [a_til zeros(1,Max-(length(a_til)-1))]
b = [b zeros(1,Max-(length(b)-1))]
c = [c zeros(1,Max-(length(c)-1))]

% State-Space cannonical observable form
A = [-a_til(2:end)' [eye(length(a_til)-2); zeros(1,length(a_til)-2)]]
B = [b(2:end)']
D = [c(2:end)'-a_til(2:end)']
H = [1 zeros(1,Max-1)]

% Prediction Matrix
for i=1:N
   G(i,:) = H*A^(i-1)*B;
end
for j=2:N
    G = [G [0; G(1:end-1,j-1)]];
end
% % for k=Nu+1:N
% %     G(:,Nu) = G(:,Nu) + G(:,k);
% % end
G = G(:,1:Nu)

for i=1:N
    F(i,:) = H*A^i;
end
F

for i=1:N
    E(i,1) = H*A^(i-1)*D;
end
E
%return
%% Simulation parameters
% Time parameters
Tsim = 60;
qntd_k = ceil(Tsim/Ts);

% Reference parameters
r1 = 0.5;
Tr1 = ceil(1/Ts);
r2 = 0.3;
Tr2 = ceil(45/Ts);

% Disturbance parameters
pert = -0.1;
Tpert = ceil(30/Ts);

% Reference and disturbance arrays
ref = [zeros(Tr1,1); r1*ones(Tr2-Tr1,1); (r1+r2)*ones(qntd_k-Tr2,1)];
q = [zeros(Tpert,1); pert*ones(qntd_k-Tpert,1)];

%% GPC
fprintf('\n GPC \n')
tic
% Quadratic equation HH
HH = 2*(G'*Q*G + lambda);
Hw = [HH 0;
    0 2*lambda_e];

% Constraints matrix Ac
Acon = [G -ones(size(G));
    -G -ones(size(G));
    1 0;
    -1 0;
    1 0;
    -1 0];

% Inicial conditions
f = zeros(N,1);
r = zeros(N,1);
y0 = 0;
du0 = 0;
u0 = 0;
x0 = zeros(length(A), 1);
e0 = 0;
[~,z0] = filter(b(2:end),a,0);

for k=1:qntd_k

    y_GPC(k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(k) = y_GPC(k) - H*x(:,k);
    
    % free response
    for i=1:N
       f(i,1) = F(i,:)*x(:,k) + E(i,1)*e(k); 
    end
    
    % reference
    r(1:N,1) = ref(k);
    
    % quadratic equation bb
    bb = 2*(f-r)'*Q*G;
    bw = [bb 0]';

    % quadratic equation c0
    c0 = (f-r)'*((f-ones(fliplr(size(r)))*r));
    
    % constraints matrix Bc
    Bcon = [y_max-f;
        -y_min+f;
        (u_max-u0);
        -(u_min-u0);
        du_max;
        -du_min];

    % QP solver
    % opt = optimoptions('quadprog','Display','off');
    % Acon = []; Bcon = []; % -> uncomment to remove constraints
    % sol = quadprog(Hw,bw,Acon,Bcon,[],[],[],[],[],opt);
    sol = quadprog(Hw,bw,Acon,Bcon);

    du_GPC(k) = sol(1);
    u_GPC(k) = u0 + du_GPC(k);

    if k>=Tpert
        [y0,z0] = filter(b(2:end),a,u_GPC(k)+pert,z0);
    else
        [y0,z0] = filter(b(2:end),a,u_GPC(k),z0);
    end

    du0 = du_GPC(k);
    u0 = u_GPC(k);
    x0 = x(:,k);
    e0 = e(k);
end
toc
% return
%% Proposed controller
fprintf('\n GPC-based PID \n')
tic
% quadratic equation HH
HH = 2*(G'*Q*G + lambda);

% constraints matrix Ac
Acon = [G -ones(size(G));
    -G -ones(size(G));
    1 0;
    -1 0;
    1 0;
    -1 0];

% initialise variables
H_T = [1 0;
    0 sqrt(HH/(2*lambda_e))];
A_til_c = Acon*H_T;

% inicial conditions
f = zeros(N,1);
r = zeros(N,1);
y0 = 0;
du0 = 0;
u0 = 0;
x0 = zeros(length(A), 1);
e0 = 0;
[~,z0] = filter(b(2:end),a,0);

for k=1:qntd_k
    % GPC
    y_PID(k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(k) = y_PID(k) - H*x(:,k);
    for i=1:N
       f(i,1) = F(i,:)*x(:,k) + E(i,1)*e(k); 
    end
    r(1:N,1) = ref(k);
    bb = 2*(f-r)'*Q*G;

    Bcon = [y_max-f;
        -y_min+f;
        (u_max-u0);
        -(u_min-u0);
        du_max;
        -du_min];

    % compute alpha_i and gama_i foi i=1 to 2N
    for i=1:2*N
        alpha(i) = -A_til_c(i,1)/A_til_c(i,2);
        beta(i) = Bcon(i)/A_til_c(i,2);
        gama(i) = 1 + alpha(i)^2;
    end

    % compute ∆u*_UC
    du_ast_UC(k) = -bb/HH;

    % compute Ug_min and Ug_max
    Ug_min = max(du_min,(u_min-u0));
    Ug_max = min(du_max,(u_max-u0));

    e_til_max = -10^10; % or -1e10

    for i=1:2*N
        e_til(i,k) = alpha(i)*du_ast_UC(k) + beta(i);
        if e_til(i,k)>e_til_max
            e_til_max = e_til(i,k);
            i_max = i;
        end
    end

    if e_til_max>0
        du_ast_SO(k) = (du_ast_UC(k) - alpha(i_max)*beta(i_max))/gama(i_max);
        du_aux = du_ast_SO(k);
    else
        du_aux = du_ast_UC(k);
    end

    du_PID(k) = du_aux;
    if du_aux<Ug_min
        du_PID(k) = Ug_min;
    end
    if du_aux>Ug_max
        du_PID(k) = Ug_max;
    end
    u_PID(k) = du_PID(k) + u0;

    % apply u(k) to the plant
    if k>=Tpert
        [y0,z0] = filter(b(2:end),a,u_PID(k)+pert,z0);
    else
        [y0,z0] = filter(b(2:end),a,u_PID(k),z0);
    end

    % update the variables
    du0 = du_PID(k);
    u0 = u_PID(k);
    x0 = x(:,k);
    e0 = e(k);
end
toc
%return
%% PID with mapping
fprintf('\n PID with mapping \n')
tic
% unconstrained controller gains
K = inv(G'*Q*G+lambda)*G'*Q;
K = K(1,:);
Kr = sum(K);
% KF = acker(A,B,[1 1 1]*0.7)
% Kr = inv(H*inv(eye(3)-A+B*KF)*B)
KF = K*F;
KE = K*E;

% augmented system
Ac = A-D*H+B*(KE*H-KF);
Bc = [B*Kr,(D-B*KE)];
Cc = KE*H-KF;
Dc = [Kr,-KE];

Controlador = tf(ss(Ac,Bc,Cc,Dc,Ts));
[nC1,dC1] = tfdata(Controlador(1),'v');
[nC2,dC2] = tfdata(Controlador(2),'v');

C_PID =  -zpk(minreal(tf(nC2,conv(dC2(1:end-1),delta),Ts,'Variable','z^-1')))
F_PID =  -zpk(minreal(tf([nC1 0],nC2,Ts,'Variable','z^-1')))

% Simulation in Simulink
out = sim('PID_Lucian2024.slx');
y_PIDmapping = out.y_PIDmapping(2:end);
u_PIDmapping = out.u_PIDmapping(2:end);
du_PIDmapping = out.du_PIDmapping(2:end);
toc
%return
%% Figure
t = Ts:Ts:Tsim;

figure
subplot(3,1,1)
hold on
stairs(t,y_GPC, 'm', 'linewidth', 2)
stairs(t,y_PIDmapping, 'b', 'linewidth', 2)
stairs(t,y_PID, '-.g', 'linewidth', 2)
stairs(t,ref, '--r', 'linewidth', 2)
stairs(t,y_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,y_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('y')
legend('GPC', 'PID with mapping', 'Proposed PID', 'Reference', 'Constraints', 'location', 'best')
title('Output')
subplot(3,1,2)
hold on
stairs(t,u_GPC, 'm', 'linewidth', 2)
stairs(t,u_PIDmapping, 'b', 'linewidth', 2)
stairs(t,u_PID, '-.g', 'linewidth', 2)
stairs(t,u_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,u_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('u')
title('Control signal')
subplot(3,1,3)
hold on
stairs(t,du_GPC, 'm', 'linewidth', 2)
stairs(t,du_PIDmapping, 'b', 'linewidth', 2)
stairs(t,du_PID, '-.g', 'linewidth', 2)
stairs(t,du_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,du_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('$\Delta u$')
title('Control signal rate')
sgtitle("\textbf{Output and Control Signals for N = }" + N + "\textbf{ and Nu = }" + Nu)

fprintf('\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')
return
% fim
%% Teste
tic
sys1 = series(C_PID, Pz);
sys2 = feedback(sys1, 1);
sys3 = series(F_PID, sys2);
[y_teste, t_out] = lsim(sys3, ref, t);
sys4 = feedback(Pz, C_PID);
[y_teste2, t_out] = lsim(sys4, q, t);
y_teste3 = y_teste + y_teste2;
t_out = t_out - Ts;

figure
hold on
stairs(t_out,y_teste3, 'b', 'linewidth', 2)
stairs(t,y_PIDmapping, '--g', 'linewidth', 2)
stairs(t,ref, '--r', 'linewidth', 2)
stairs(t,y_max*ones(qntd_k,1), '-.k', 'linewidth', 1)
stairs(t,y_min*ones(qntd_k,1), '-.k', 'linewidth', 1)
% grid on
% axis tight
xlabel('Time (s)')
ylabel('y')
legend('without Simulink','with Simulink', 'Reference', 'Constraints', 'location', 'best')
title('Output')
toc

% %
% [r_teste, t_out] = lsim(F_PID, ref, t);
% e_teste = r_teste - y_teste3;
% [u_teste, t_out] = lsim(C_PID, e_teste, t);
% t_out = t_out - Ts;
% figure
% stairs(t_out,u_teste, 'r', 'linewidth', 2)
% %