function Blanca2023_solver(block)

setup(block);

function setup(block)

% Register number of dialog parameters
block.NumDialogPrms = 16;
N = block.DialogPrm(1).Data;
Nu = block.DialogPrm(2).Data;
P = block.DialogPrm(3).Data;
V_B = block.DialogPrm(4).Data;
G_eq = block.DialogPrm(5).Data;
P_eq = block.DialogPrm(6).Data;
V_B_eq = block.DialogPrm(7).Data;
u_max = block.DialogPrm(8).Data;
u_min = block.DialogPrm(9).Data;
du_max = block.DialogPrm(10).Data;
du_min = block.DialogPrm(11).Data;
H_aug = block.DialogPrm(12).Data;
n_eq = block.DialogPrm(13).Data;
Ts = block.DialogPrm(14).Data;
G = block.DialogPrm(15).Data;
gamma = block.DialogPrm(16).Data;

% [ny,~] = size(y_max);
[~,nu] = size(u_max);
% [~,ma] = size(A);

% Register number of ports
block.NumInputPorts  = 4;
block.NumOutputPorts = 2;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Dimensions        = 2;
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

block.InputPort(4).Dimensions        = nu;
block.InputPort(4).DatatypeID  = 0;  % double
block.InputPort(4).Complexity  = 'Real';
block.InputPort(4).DirectFeedthrough = true;

% Override output port properties
block.OutputPort(1).Dimensions       = nu;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';

block.OutputPort(2).Dimensions       = nu;
block.OutputPort(2).DatatypeID  = 0; % double
block.OutputPort(2).Complexity  = 'Real';

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [Ts 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Terminate', @Terminate); % Required

function Outputs(block)
N = block.DialogPrm(1).Data;
Nu = block.DialogPrm(2).Data;
P = block.DialogPrm(3).Data;
V_B = block.DialogPrm(4).Data;
G_eq = block.DialogPrm(5).Data;
P_eq = block.DialogPrm(6).Data;
V_B_eq = block.DialogPrm(7).Data;
u_max = block.DialogPrm(8).Data;
u_min = block.DialogPrm(9).Data;
du_max = block.DialogPrm(10).Data;
du_min = block.DialogPrm(11).Data;
H_aug = block.DialogPrm(12).Data;
n_eq = block.DialogPrm(13).Data;
Ts = block.DialogPrm(14).Data;
G = block.DialogPrm(15).Data;
gamma = block.DialogPrm(16).Data;

% ------------------ main ----------------------
I_L = tril(ones(Nu, Nu));

A_du = [eye(Nu); -eye(Nu)];
b_du = [du_max * ones(Nu, 1); -du_min * ones(Nu, 1)];
A_u = [I_L; -I_L];
A_ineq = [A_du; A_u];
A_ineq_aug = [A_ineq, zeros(size(A_ineq, 1), n_eq)];


% Input parameters
r=ones(N,1)*block.InputPort(1).Data;
u_r = block.InputPort(2).Data;
x_m = block.InputPort(3).Data;
u0_prev = block.InputPort(4).Data;

f = P * x_m + V_B * u0_prev;

f_du = 2 * G' * gamma * (f - r);

% Augmented linear objective function f (zeros for slack variables cost)
f_aug = [f_du; zeros(n_eq, 1)];

% G_eq * Delta_u0 - xi_2 = y_ref_eq - (P_eq * x_m + V_B_eq * u0_prev)
Aeq_aug = [G_eq, -eye(n_eq)]; 
beq_aug = r(1) * ones(n_eq, 1) - (P_eq * x_m + V_B_eq * u0_prev);
    
b_u_upper =  u_max * ones(Nu, 1) - u0_prev * ones(Nu, 1) + u_r * ones(Nu, 1);
b_u_lower = -u_min * ones(Nu, 1) + u0_prev * ones(Nu, 1) - u_r * ones(Nu, 1);
b_u = [b_u_upper; b_u_lower];
b_ineq = [b_du; b_u];

Delta_U_aug = quadprog(H_aug, f_aug, A_ineq_aug, b_ineq, Aeq_aug, beq_aug, [], [], [], []);

delta_u0_k = Delta_U_aug(1);
u0_k = u0_prev + delta_u0_k;

block.OutputPort(1).Data = u0_k;
block.OutputPort(2).Data = delta_u0_k;


function Terminate(block)


