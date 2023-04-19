% Parameters
fs = 10e3; % Sampling frequency 10kHz
Ts = 1/fs; % Sampling time 100ns
omega_sync = 314.159; % synchronous speed 314.159rad/s (3000 rpm)
Ld = 0.8e-3; % d-axis inductance 0.8mH
Lq = 0.8e-3; % q-axis inductance 0.8mH
Rs = 45e-3; % stator resistance 45mΩ
psi_ro = 0.65; % Magnetic flux linkage 0.65Wb
Ir = 0.35; % Moment of inertia 0.35 kg.m^2
Br = 0.004; % Viscosity coefficient

Q = [100 0 0; 0 1 0; 0 0 30];
R = [0.8 0; 0 0.8];

tl_mean = 0; % load torque mean

omegass = omega_sync;
idss = 0;
iqss = 2*(Br*omegass + tl_mean)/(3*psi_ro);
vdss = - Lq * omegass * iqss;
vqss = omegass * psi_ro + Rs * iqss;

A = [(1 - Ts*Rs/Ld), (Ts*Lq*omegass/Ld), Ts*iqss; -(Ts*Ld*omegass/Lq), (1 - Ts*Rs/Lq), -Ts*(psi_ro/Lq + idss); 0, 1.5*Ts*psi_ro/Ir, (1 - Ts*Br/Ir)];
B = [Ts/Ld, 0; 0, Ts/Lq; 0, 0];

wbar = [0; 0; 0]; % noise mean

% for Gaussian with stddev = 1
% tl_sd = 1;
% M3 = [0; 0; 0];
% m4 = 1800;
% for exponential with gamma = 1 (shifted to be zero-mean)
% tl_sd = 1;
% M3 = [0; 0; 60];
% m4 = 7200;
% for uniform from -1 to 1
tl_sd = sqrt(1/3);
M3 = [0; 0; 0];
m4 = 80;
% for pareto with alpha = 5 (shifted to be zero-mean)
% tl_sd = sqrt(5/48);
% M3 = [0; 0; 4.6875];
% m4 = 710.9375;

W = [0, 0, 0; 0, 0, 0; 0, 0, tl_sd^2]; % noise covariance

num_states = 1000; % control sequence length

% Gaussian: -3, 1
mu_array = logspace(-6,6,100);
epsilon_array = zeros(1,100);
for mu_count = 1:100
    mu = mu_array(mu_count);
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
        Vt = V(:,:,n+1);
        St = S(:,:,n+1);
        Tt = T(:,:,n+1);
        temp = inv(B'*Vt*B+R);
        V(:,:,n) = A'*Vt*A + Q_mu - A'*Vt*B*temp*B'*Vt*A;
        K(:,:,n) = -temp*B'*Vt*A;
        S(:,:,n) = (A+B*K(:,:,n))'*St+Q;
        T(:,:,n) = (A+B*K(:,:,n))'*(Tt+Vt);
        l(:,n) = -2*mu*temp*B'*St*M3;
        h(:,n) = -temp*B'*(Tt+Vt)*wbar;
    end

    P = zeros(3,3,num_states);
    z = zeros(3,num_states);
    r = zeros(1,num_states);

    P(:,:,num_states) = 4*Q*W*Q;
    z(:,num_states) = 4*M3'*Q;

    for n = num_states-1:-1:1
        Pt = P(:,:,n+1);
        zt = z(:,n+1);
        temp1 = A + B*K(:,:,n);
        temp2 = (B*l(:,n) + B*h(:,n) + wbar);
        P(:,:,n) = temp1'*Pt*temp1 + 4*Q*W*Q;
        z(:,n) = temp1'*zt + 4*Q*M3 + 2*temp1'*Pt*temp2;
        r(n) = r(n+1) + trace(Pt*W) + zt'*temp2 + temp2'*Pt*temp2;
    end

    w = randn(1,num_states) * tl_sd;
    x = zeros(3,num_states);
    u = zeros(2,num_states);

    epsilon_bar = x(:,1)'*P(:,:,1)*x(:,1) + z(:,1)'*x(:,1) + r(1);
    epsilon = epsilon_bar + num_states*m4 - 4*num_states*trace(Q'*W'*W*Q);
    epsilon_array(mu_count) = epsilon;

    for n = 1:num_states-1
        u(:,n) = K(:,:,n)*x(:,n) + l(:,n) + h(:,n);
        x(:,n+1) = A*x(:,n) + B*u(:,n) + w(:,n);
    end
end
figure
semilogx(mu_array,epsilon_array)