function PID_GPC(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.

%   Copyright 2003-2018 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C MEX counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of dialog parameters
block.NumDialogPrms = 16;
F = block.DialogPrm(1).Data;
E = block.DialogPrm(2).Data;
G = block.DialogPrm(3).Data;
y_max = block.DialogPrm(4).Data;
y_min = block.DialogPrm(5).Data;
u_max = block.DialogPrm(6).Data;
u_min = block.DialogPrm(7).Data;
du_max = block.DialogPrm(8).Data;
du_min = block.DialogPrm(9).Data;
Ts = block.DialogPrm(10).Data;
N = block.DialogPrm(11).Data;
Nu = block.DialogPrm(12).Data;
Q = block.DialogPrm(13).Data;
lambda = block.DialogPrm(14).Data;
lambda_e = block.DialogPrm(15).Data;
A = block.DialogPrm(16).Data;

[ny,~] = size(y_max);
[~,nu] = size(u_max);
[~,ma] = size(A);
% Register number of ports
block.NumInputPorts  = 4;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = ma;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions        = ny;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Dimensions        = nu;
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

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
%%
function Outputs(block)
F = block.DialogPrm(1).Data;
E = block.DialogPrm(2).Data;
G = block.DialogPrm(3).Data;
y_max = block.DialogPrm(4).Data;
y_min = block.DialogPrm(5).Data;
u_max = block.DialogPrm(6).Data;
u_min = block.DialogPrm(7).Data;
du_max = block.DialogPrm(8).Data;
du_min = block.DialogPrm(9).Data;
Ts = block.DialogPrm(10).Data;
N = block.DialogPrm(11).Data;
Nu = block.DialogPrm(12).Data;
Q = block.DialogPrm(13).Data;
lambda = block.DialogPrm(14).Data;
lambda_e = block.DialogPrm(15).Data;
A = block.DialogPrm(16).Data;


% T = tril(ones(Nu));
% b0 = [1;zeros(Nu-1,1)];

% Input parameters
r=ones(N,1)*block.InputPort(4).Data;
u0 = block.InputPort(3).Data;
x = block.InputPort(1).Data;
e = block.InputPort(2).Data;

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

f=zeros(N,1);
% free response
for i=1:N
    f(i,1) = F(i,:)*x + E(i,1)*e; 
end

% quadratic equation bb
bb = 2*(f-r)'*Q*G;

% constraints matrix Bc
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
du_ast_UC = -bb/HH;

% compute Ug_min and Ug_max
Ug_min = max(du_min,(u_min-u0));
Ug_max = min(du_max,(u_max-u0));

e_til_max = -10^10; % or -1e10

for i=1:2*N
    e_til(i) = alpha(i)*du_ast_UC + beta(i);
    if e_til(i)>e_til_max
        e_til_max = e_til(i);
        i_max = i;
    end
end

if e_til_max>0
    du_ast_SO = (du_ast_UC - alpha(i_max)*beta(i_max))/gama(i_max);
    du_aux = du_ast_SO;
else
    du_aux = du_ast_UC;
end

% saturation of Δu
du = du_aux;
if du_aux<Ug_min
    du = Ug_min;
end
if du_aux>Ug_max
    du = Ug_max;
end

block.OutputPort(1).Data = du;

%end Outputs

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate
