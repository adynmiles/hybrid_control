%% Hybrid Optimization Plot %%

num_impulses = 0:5:40;
num_impulses_interp = 0:1:40;
settling_time = [23, 15, 9.5, 8, 6, 4.5, 3.5, 3, 3];  
input_current = [3.19, 3.21, 3.21, 3.47, 3.21, 3.53, 3.93, 4.09, 3.94];
mag_torque_norm = [0.0073, 0.0095, 0.0072, 0.0065, 0.0041, 0.0032, 0.0055, 0.0025, 0.0024];
imp_torque_norm = [0, 0.00024184, 0.00031629, 0.00039259, 0.00039691, 0.00046640, 0.00066626, 0.00060925, 0.00055208];
fig_of_merit = settling_time.*(mag_torque_norm.^2 + imp_torque_norm.^2);
fig_of_merit_interp = 1000*interp1(num_impulses, fig_of_merit, num_impulses_interp, 'linear');       % Interpolate available power



figure;
plot(num_impulses_interp, fig_of_merit_interp, 'LineWidth', 2);
title('Number of Impulses Optimization');
ylabel('Figure of Merit');
xlabel('Number of Impulses over 10 orbits');

