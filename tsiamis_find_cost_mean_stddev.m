% clear all
% 
% load('data.mat')
% 
% mu = 10^(0.5);
% 
% id = ans(2,:);
% iq = ans(3,:);
% omega = ans(4,:);
% vd = ans(5,:);
% vq = ans(6,:);
% 
% save('tsiamis_05_simulink.mat')
%% 
% len = length(id);
% 
% cost = zeros(1,len);
% omega_cost = cost;
% x = [id; iq; omega];
% u = [vd; vq];
% Q = diag([100, 1, 30]);
% R = diag([0.8,0.8]);
% for count = 1:len
%     cost(count) = x(:,count)'*Q*x(:,count) + u(:,count)'*R*u(:,count);
%     omega_cost(count) = omega(count)^2;
% end
% 
% curr_mean = sum(cost)/len
% curr_stddev = sqrt(sum(cost.^2)/len - curr_mean^2)
% curr_omega_mean = sum(omega_cost)/len
% curr_omega_stddev = sqrt(sum(omega_cost.^2)/len - curr_omega_mean^2)
%%
% load('simulink_tsiamis_lqr.mat')
% 
% mu_array(end+1) = mu;
% tsiamis_mean(end+1) = curr_mean;
% tsiamis_stddev(end+1) = curr_stddev;
% omega_mean(end+1) = curr_omega_mean;
% omega_stddev(end+1) = curr_omega_stddev;
% 
% save('simulink_tsiamis_lqr.mat','mu_array', 'tsiamis_stddev', 'tsiamis_mean','omega_stddev','omega_mean')

%% 
load('simulink_tsiamis_lqr.mat')
figure ('position',[100 100 1600 300])
subplot(1,3,2)
semilogx(mu_array(1,2:end), tsiamis_mean(1,2:end), 'bo--')
hold on
% plot(mu_array(1,1), tsiamis_mean(1,1), 'r*')
grid on
% xlim([-0.5 4])

subplot(1,3,1)
semilogx(mu_array(1,2:end), omega_stddev(1,2:end), 'bo--')
hold on
% plot(mu_array(1,1), omega_stddev(1,1), 'r*')
grid on
xlim([-0.5 4])

subplot(1,3,3)
plot(tsiamis_mean(1,2:end), omega_stddev(1,2:end), 'bo--')
hold on
% plot(tsiamis_mean(1,1), omega_stddev(1,1), 'r*')
grid on