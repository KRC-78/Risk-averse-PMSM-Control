clear
close all
clc
%%
% ABB HDS180-4876B Parameters
f = 1e3; % Sampling frequency 1kHz
T = 1/f; % Sampling time 100ns
omega_rpm = 1500;
omega_sync = omega_rpm*2*pi/60; % synchronous speed
Ld = 1.45e-3; % d-axis inductance
Lq = 1.45e-3; % q-axis inductance
Rs = 65e-3; % stator resistance
psi_ro = 0.16513; % Magnetic flux linkage
Ir = 0.00821; % Moment of inertia
Br = 0.004; % Viscosity coefficient
p = 5; % # of pole pairs

Q = [100 0 0; 0 1 0; 0 0 30];
R = [0.8 0; 0 0.8];

tl_mean = 0; % load torque mean
tl_sd = 1; % load torque standard deviation
tl_var = tl_sd^2;

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*p*psi_ro);
vdss = - Lq * omegass * p*iqss;
vqss = omegass * p*psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*p*Lq*omegass/Ld), T*p*iqss; -(T*Ld*p*omegass/Lq), (1 - T*Rs/Lq), -T*p*(psi_ro/Lq + idss); 0, 1.5*T*p*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];
D = [0; 0; -T/Ir];

wbar = [0; 0; 0]; % noise mean
W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];

tl_sd = 1; % load torque standard deviation
M3 = [0; 0; 0];
m4 = 1800;
%% Tsiamis

mu = 10^(0.5);
% min mu 1e-1.5 max mu 1e0.5 incremental 0.25, no Ts/Ir coefficient in W
Q_mu = Q + 4*mu*Q*W*Q;

V = zeros(3,3,num_states);
K = zeros(2,3,num_states);
S = zeros(3,3,num_states);
T = zeros(3,3,num_states);
l = zeros(2,num_states);
h = zeros(2,num_states);

V(:,:,num_states) = Q_mu;
S(:,:,num_states) = Q;

for n = num_states-1:-1:1
    V_prev = V(:,:,n+1);
    S_prev = S(:,:,n+1);
    T_prev = T(:,:,n+1);
    temp = inv(B'*V_prev*B+R);
    V(:,:,n) = A'*V_prev*A + Q_mu - A'*V_prev*B*temp*B'*V_prev*A;
    K(:,:,n) = -temp*B'*V_prev*A;
    S(:,:,n) = (A+B*K(:,:,n))'*S_prev+Q;
    T(:,:,n) = (A+B*K(:,:,n))'*(T_prev+V_prev);
    l(:,n) = -2*mu*temp*B'*S_prev*M3;
    h(:,n) = -temp*B'*(T_prev+V_prev)*wbar;
end

controller_gain = K(:,:,1);
controller_offset = l(:,1) + h(:,1);
