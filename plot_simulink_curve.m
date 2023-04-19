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

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*p*psi_ro);
vdss = - Lq * omegass * p*iqss;
vqss = omegass * p*psi_ro + Rs * iqss;

A = [(1 - T*Rs/Ld), (T*p*Lq*omegass/Ld), T*p*iqss; -(T*Ld*p*omegass/Lq), (1 - T*Rs/Lq), -T*p*(psi_ro/Lq + idss); 0, 1.5*T*p*psi_ro/Ir, (1 - T*Br/Ir)];
B = [T/Ld, 0; 0, T/Lq; 0, 0];

wbar = [0; 0; 0]; % noise mean
W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 2000; % control sequence length, 0.2 secs

xref = zeros(3,num_states);
uref = zeros(2,num_states);

xss = [idss; iqss; omegass];
uss = [vdss; vqss];

load('lqr_simulink.mat')

t = linspace(0, 10, length(id));

figure ('position',[100 100 1600 300])
subplot(1,3,1)
plot(t, id)
hold on

subplot(1,3,2)
plot(t, iq + iqss)
hold on

subplot(1,3,3)
plot(t, omega + omegass)
hold on

load('leqg_045_simulink.mat')

t = linspace(0, 10, length(id));

subplot(1,3,1)
plot(t, id)
grid on
legend('LQR','LEQG')
% legend('Standard', 'Constrained')
ylim([-3.3e-2 3.3e-2])

subplot(1,3,2)
plot(t, iq + iqss)
grid on
legend('LQR','LEQG')
% legend('Standard', 'Constrained')
ylim([-1 2.5])

subplot(1,3,3)
plot(t, omega + omegass)
grid on
legend('LQR','LEQG')
% legend('Standard', 'Constrained')
ylim([156.3 158])