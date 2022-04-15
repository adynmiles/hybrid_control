%% Hybrid Continuous Comparison %%
load continuous_1.mat
load continuous_time_1.mat
load continuous_current_1.mat
t_continuous = t_out;
theta_continuous = theta_out;
current_continuous = I;
load hybrid_1a.mat
load hybrid_time_1a.mat
load hybrid_current_1a.mat
t_hybrid1a = t_full;
theta_hybrid1a = theta_full;
current_hybrid1a = I;
load hybrid_1b.mat
load hybrid_time_1b.mat
load hybrid_current_1b.mat
t_hybrid1b = t_full;
theta_hybrid1b = theta_full;
current_hybrid1b = I;
load hybrid_1c.mat
load hybrid_time_1c.mat
load hybrid_current_1c.mat
t_hybrid1c = t_full;
theta_hybrid1c = theta_full;
current_hybrid1c = I;


figure;
subplot(3, 1, 1); 
plot(t_continuous/5400, theta_continuous(:,1), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, theta_hybrid1a(:,1), 'LineWidth', 2);
plot(t_hybrid1b/5400, theta_hybrid1b(:,1), 'LineWidth', 2);
plot(t_hybrid1c/5400, theta_hybrid1c(:,1), 'LineWidth', 2);
title("Roll Behaviour Comparison for Hybrid and Continuous Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);
legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)

subplot(3, 1, 2); 
plot(t_continuous/5400, theta_continuous(:,2), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, theta_hybrid1a(:,2), 'LineWidth', 2);
plot(t_hybrid1b/5400, theta_hybrid1b(:,2), 'LineWidth', 2);
plot(t_hybrid1c/5400, theta_hybrid1c(:,2), 'LineWidth', 2);
title("Pitch Behaviour Comparison for Hybrid and Continuous Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);
%legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)

subplot(3, 1, 3); 
plot(t_continuous/5400, theta_continuous(:,3), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, theta_hybrid1a(:,3), 'LineWidth', 2);
plot(t_hybrid1b/5400, theta_hybrid1b(:,3), 'LineWidth', 2);
plot(t_hybrid1c/5400, theta_hybrid1c(:,3), 'LineWidth', 2);
title("Yaw Behaviour Comparison for Hybrid and Continuous Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);
%legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)


figure;
subplot(3, 1, 1); 
plot(t_continuous/5400, current_continuous(:,1), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, current_hybrid1a(:,1), 'LineWidth', 2);
plot(t_hybrid1b/5400, current_hybrid1b(:,1), 'LineWidth', 2);
plot(t_hybrid1c/5400, current_hybrid1c(:,1), 'LineWidth', 2);
title("Control Input Current Required (Roll)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current (A)", 'FontSize', 14);
%legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)

subplot(3, 1, 2); 
plot(t_continuous/5400, current_continuous(:,2), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, current_hybrid1a(:,2), 'LineWidth', 2);
plot(t_hybrid1b/5400, current_hybrid1b(:,2), 'LineWidth', 2);
plot(t_hybrid1c/5400, current_hybrid1c(:,2), 'LineWidth', 2);
title("Control Input Current Required (Pitch)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current (A)", 'FontSize', 14);
legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)

subplot(3, 1, 3); 
plot(t_continuous/5400, current_continuous(:,3), 'LineWidth', 2);
hold on;
plot(t_hybrid1a/5400, current_hybrid1a(:,3), 'LineWidth', 2);
plot(t_hybrid1b/5400, current_hybrid1b(:,3), 'LineWidth', 2);
plot(t_hybrid1c/5400, current_hybrid1c(:,3), 'LineWidth', 2);
title("Control Input Current Required (Yaw)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current (A)", 'FontSize', 14);
%legend('Continuous', 'Hybrid - 10 Impulses', 'Hybrid - 20 Impulses', 'Hybrid - 30 Impulses', 'FontSize', 14)