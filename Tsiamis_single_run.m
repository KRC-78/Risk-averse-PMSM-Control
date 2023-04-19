clear
close all
%%

distribution_type = "gaussian";

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

% All noise distributions are zero-mean
wbar = [0; 0; 0]; % noise mean

if distribution_type == "gaussian" % for Gaussian with stddev = 1
    tl_sd = 1; % load torque standard deviation
    M3 = [0; 0; 0];
    m4 = 1800;
    pd = makedist('Normal','mu',0,'sigma',1);
elseif distribution_type == "logistic" % for logistic with mu = 1
    tl_sd = pi / sqrt(3);
    M3 = [0; 0; 0];
    m4 = 320 * pi^4;
    pd = makedist('Logistic','mu',0,'sigma',1);
elseif distribution_type == "uniform" % for uniform from -1 to 1
    tl_sd = sqrt(1/3);
    M3 = [0; 0; 0];
    m4 = 80;
    pd = makedist('Uniform','lower',-1,'upper',1);
elseif distribution_type == "tstudent" % for t-student with nu = 5
    tl_sd = sqrt(5/3);
    M3 = [0; 0; 0];
    m4 = 900*(gamma(2.5)*gamma(0.5)*25/(sqrt(pi)*gamma(2.5)) - gamma(1.5)*gamma(1.5)*5/(sqrt(pi)*gamma(2.5)));
    pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',5);
end

W = [0, 0, 0; 0, 0, 0; 0, 0, (Ts/Ir)^2 * tl_sd^2]; 
% noise covariance, thesis version lacks Ts/Ir coefficient
% however wouldn't change qualitative result, only requires larger mu values

num_states = 50000; % control sequence length

xss = [idss; iqss; omegass];
uss = [vdss; vqss];
%% Risk constrained LQR
mu_nums = 20;
min_mu = 5; % Gaussian -3, logistic -3, uniform -2, t-student -2
max_mu = 10; % Gaussian 1, logistic 1, uniform 1, t-student 2
mu_array = logspace(min_mu,max_mu,mu_nums);

omega_cost = zeros(num_states, mu_nums);
omega_cost_mean = zeros(1,mu_nums);
omega_cost_stddev = zeros(1,mu_nums);
x_cost = zeros(num_states, mu_nums);
x_cost_mean = zeros(1,mu_nums);
x_cost_stddev = zeros(1,mu_nums);

for mu_count = 1:mu_nums
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

    x = zeros(3,num_states);
    u = zeros(2,num_states);
    tl_noise = random(pd,1,num_states);

    for n = 1:num_states-1
        u(:,n) = K(:,:,n) * x(:,n) + l(:,n) + h(:,n);
        x(:,n+1) = A*x(:,n) + B*u(:,n) + [0;0;-Ts*tl_noise(n)/Ir];
        omega_cost(n, mu_count) = x(3,n)^2;
        x_cost(n,mu_count) = x(:,n)'*Q*x(:,n);
    end
    omega_cost(num_states,mu_count) = x(3, num_states)^2;
    omega_cost_mean(mu_count) = sum(omega_cost(:,mu_count))/num_states;
    omega_cost_stddev(mu_count) = sqrt(sum(omega_cost(:,mu_count).^2)/num_states - omega_cost_mean(mu_count)^2);
    x_cost(num_states,mu_count) = x(:,num_states)'*Q*x(:,num_states);
    x_cost_mean(mu_count) = sum(x_cost(:,mu_count))/num_states;
    x_cost_stddev(mu_count) = sqrt(sum(x_cost(:,mu_count).^2)/num_states - x_cost_mean(mu_count)^2);
end

%% Plot state cost mean vs omega stddev
figure
hold on
grid on
title('State cost mean vs omega cost stddev')

plot(x_cost_mean, omega_cost_stddev)

xlabel('State mean')
ylabel('Omega stddev')

%% 
% if distribution_type == "gaussian" % for Gaussian with stddev = 1
%     save('tsiamis_gaussian.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
% elseif distribution_type == "logistic" % for logistic with mu = 1
%     save('tsiamis_logistic.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
% elseif distribution_type == "uniform" % for uniform from -1 to 1
%     save('tsiamis_uniform.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
% elseif distribution_type == "tstudent" % for t-student with nu = 5
%     save('tsiamis_tstudent.mat', "x_cost_mean", "omega_cost_stddev", "mu_array", "num_states")
% end


