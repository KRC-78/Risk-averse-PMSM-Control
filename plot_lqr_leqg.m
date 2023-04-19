tl_sd = 2

if tl_sd == 0.5
    load('lqr_gaussian_05.mat')
    load('leqg_gaussian_05.mat')
    theta_array = linspace(-0.3, -3, 10);
    min_theta = -3.6;
    max_theta = 0.6;
elseif tl_sd == 1
    load('lqr_gaussian_1.mat')
    load('leqg_gaussian_1.mat')
    theta_array = linspace(-0.1, -1, 10);
    min_theta = -1.2;
    max_theta = 0.2;
elseif tl_sd == 2
    load('lqr_gaussian_2.mat')
    load('leqg_gaussian_2.mat')
    theta_array = linspace(-0.02, -0.22, 10);
    min_theta = -0.26;
    max_theta = 0.04;
end
lqr_curr_mean = lqr_mean(end);
lqr_curr_stddev = lqr_stddev(end);

leqg_curr_mean = leqg_mean(end,:);
leqg_curr_stddev = leqg_stddev(end,:);
%% 
figure ('position',[100 100 1600 300])
subplot(1,3,1)
plot(theta_array, leqg_curr_mean, 'bo--')
hold on
plot(0, lqr_curr_mean, 'r*')
xlim([min_theta max_theta])
grid on

subplot(1,3,2)
plot(theta_array, leqg_curr_stddev, 'bo--')
hold on
plot(0, lqr_curr_stddev, 'r*')
xlim([min_theta max_theta])
grid on

subplot(1,3,3)
plot(leqg_curr_mean, leqg_curr_stddev, 'bo--')
hold on
plot(lqr_curr_mean, lqr_curr_stddev, 'r*')
grid on