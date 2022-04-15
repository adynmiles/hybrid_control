%% MAGNETIC AND IMPULSIVE SATELLITE CONTROL  %%

initial_values = deg2rad([5; 5; 5; 1; 1; 1]);               % Input initial conditions in degrees.

impulses = 30;                                               % Number of impulses.

n = 400;                                                    % Number of turns
d = 0.100;                                                  % Diameter of magnetorquer (m)
A = (pi*d^2)/4;                                             % Area of magnetorquers (m^2)     
omega = 100;
dt = 1;                                                     % Time step
T = 54000;                                                  % Total simulation time
t_list = 0:dt:T;                                            % Time grid for plotting                         
t_span = [0, T];
tk = linspace(0, T, impulses);                              % Discrete time array

inc = deg2rad(86);                                          % Inclination

G = 6.674e-11;                                              % Universal gravitational constant
M_e = 5.97e24;                                              % Mass of Earth
mu_e = G * M_e;                                             % Gravitational constant of Earth

I1 = 22;                                                    % Inertia values                           
I2 = 30;
I3 = 15;

Q_matrix = eye(6);                                          % LQR Q matrix input
Qd_matrix = eye(6);                                         % LQR Q matrix discrete input
rho = 0.1;                                                    % LQR matrix scaling factor
rho_d = 0.001;                                               % Discrete LQR scaling factor
R_matrix = rho*eye(3);                                      % LQR R matrix input
Rd_matrix = rho_d*eye(3);                                   % LQR R matrix discrete input

R_0 = 6871000;                                              % Orbital radius (assuming 500km altitude)
omega_o = sqrt(mu_e/R_0^3);                                 % Natural frequency
k_1 = (I2 - I3)/I1; 
k_3 = (I2 - I1)/I3;

I_inv = [-1/I1, 0, 0;
        0, -1/I2, 0;
        0, 0, -1/I3];

Bd = [0, 0, 0;
      0, 0, 0; 
      0, 0, 0;
      I_inv(1,1), I_inv(1,2), I_inv(1,3);
      I_inv(2,1), I_inv(2,2), I_inv(2,3);
      I_inv(3,1), I_inv(3,2), I_inv(3,3)];

Ac = zeros(6);                                              % Dynamics matrix
    Ac(1, 4) = 1;
    Ac(2, 5) = 1;
    Ac(3, 6) = 1;
    Ac(4, 1) = -4*k_1*omega_o^2;
    Ac(4, 6) = ((1 - k_1)*omega_o);
    Ac(5, 2) = -(3*(omega_o^2)*(I1 - I3));
    Ac(6, 3) = -(k_3*omega_o^2); 
    Ac(6, 4) = -((1 - k_3)*omega_o);
    
options = odeset('RelTol',1e-8,'AbsTol',1e-8);              % Error tolerance for numerical integration

% Get open loop solution to find Bc(t) at each time step.
[t_open, theta_open] = ode45(@(t_1, theta_1) open_loop(t_1, theta_1, Ac) , t_list, initial_values, options);
[Bc, b_o] = solve_Bc(t_open, theta_open);

% Initialize matrices for Riccati equation
K = zeros(length(t_list), 3, 6);
Kd = zeros(length(t_list), 3, 6);
P = zeros(length(t_list), 6, 6);
Pd = zeros(length(t_list), 6, 6);
P(end, :, :) = eye(6);                                      % Terminal condition for P
Pd(end, :, :) = eye(6);
Bc_final = squeeze(Bc(end, :, :));
P_final = squeeze(P(end, :, :));
K(end, :, :) = inv(R_matrix)*Bc_final.'*P_final;
Kd(end, :, :) = inv(R_matrix)*Bd.'*P_final; 

% Apply Euler scheme to solve Riccati equation backwards
for i = length(t_list):-1:2
    P_i = squeeze(P(i, :, :));
    Bc_i = squeeze(Bc(i, :, :));
    P_iminusone = P_i + dt*(P_i*Ac + (Ac.'*P_i) - (P_i*Bc_i*inv(R_matrix)*Bc_i.'*P_i) + Q_matrix);
    P(i - 1, :, :) = P_iminusone;
    Bc_iminusone = squeeze(Bc(i - 1, :, :));
    K(i - 1, :, :) = inv(R_matrix)*Bc_iminusone.'*P_iminusone;
end

t_full = [];
theta_full = [];

for i = impulses:-1:2
    Pd_i = squeeze(Pd(i, :, :));
    P_diminusone = Qd_matrix + Pd_i - (Pd_i*Bd*inv(Rd_matrix + (Bd.'*Pd_i*Bd))*Bd.'*Pd_i);
    Pd(i - 1, :, :) = P_diminusone;
    Kd(i - 1, :, :) = inv(Rd_matrix)*Bd.'*inv(eye(6).')*(P_diminusone - Qd_matrix);
end

theta_impulses = [];
% Using optimal gains, numerically integrate
for i = 1:impulses-1
    [t_out, theta_out] = ode45(@(t, theta) lqr_control(t, theta, Ac, Bc, K, t_list, i, T, impulses), [tk(i), tk(i+1)], initial_values, options);
    t_full = [t_full; t_out];
    theta_full = [theta_full; theta_out];
    Kd_current = squeeze(Kd(i, :, :));
    theta_impulse = (eye(6) - (Bd*Kd_current))*theta_full(end,:,:)';
    theta_impulses = [theta_impulses, theta_impulse];
    theta_full(end, :) = theta_impulse;
    initial_values = theta_impulse;
end 

% Find control inputs
m = zeros(length(t_full), 3);

for i = 1:length(t_full)
    %t_diff = abs(t_full(i) - t_list);
    %[val, idx] = min(t_diff);
    m(i, :) = squeeze(K(i, :, :))*theta_full(i, :)';
end

I = m / (n*A);                                              % Input current   
E = (omega)/((n^2) * (A^2))* (T*((m(:,1)'*m(:,1)) + (m(:,2)'*m(:,2)) + (m(:,3)'*m(:,3))));
fprintf("Total Energy Consumed: %dJ\n", E);

T = zeros(impulses, 3);
for i = 1:impulses-1
    T(i, :) = squeeze(Kd(i, :, :))*theta_impulses(:, i);
end

% Assume 60s to 100s
% Integrate the differential equation.
% 2nd order differential equation, I*theta(theta_doubledot) = delta_v
torque_imp = sqrt(impulses/54000*((T(:,1)'*T(:,1) + T(:,2)'*T(:,2) + T(:,3)'*T(:,3)))/impulses);

torque_mag = zeros(length(t_out), 3);
for i = 1:length(t_out)
    m_cross_i = [0, -m(i,3), m(i,2);
                 m(i,3), 0, -m(i,1);
                 -m(i,2), m(i,1), 0];
    t_diff = abs(t_out(i) - t_list);
    [val, idx] = min(t_diff);
    torque_mag(i, :) = m_cross_i*b_o(idx, :)';
end
    
 
total_mag = sqrt(54000*((torque_mag(:,1)'*torque_mag(:,1) + torque_mag(:,2)'*torque_mag(:,2) + torque_mag(:,3)'*torque_mag(:,3)))/(54000));
total_torque = sqrt(total_mag^2 + torque_imp^2);

% Plot angular position
figure;
subplot(3, 1, 1);
plot(t_full/5400, theta_full(:,1), 'LineWidth', 2);
title("Roll Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 2);
plot(t_full/5400, theta_full(:,2), 'LineWidth', 2);
title("Pitch Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 3);
plot(t_full/5400, theta_full(:,3), 'LineWidth', 2);
title("Yaw Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

% Plot angular velocities
figure;
subplot(3, 1, 1);
plot(t_full/5400, theta_full(:,4), 'LineWidth', 2);
title("Roll Angular Velocity Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

subplot(3, 1, 2);
plot(t_full/5400, theta_full(:,5), 'LineWidth', 2);
title("Pitch Angular Velocity Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

subplot(3, 1, 3);
plot(t_full/5400, theta_full(:,6), 'LineWidth', 2);
title("Yaw Angular Velocity Behaviour with Hybrid LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

% Plot required control input current
figure;
subplot(3, 1, 1); 
plot(t_full/5400, I(:, 1), 'LineWidth', 2);
title("Control Input Current Required (Roll)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);

subplot(3, 1, 2); 
plot(t_full/5400, I(:, 2), 'LineWidth', 2);
title("Control Input Current Required (Pitch)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);

subplot(3, 1, 3); 
plot(t_full/5400, I(:, 3), 'LineWidth', 2);
title("Control Input Current Required (Yaw)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);



%% State Space Equation with Magnetic

function [thetadot] = open_loop(t, theta, Ac)
    G = 6.674e-11;                      
    M_e = 5.97e24;                      
    mu_e = G * M_e;                     

    I1 = 22;                            
    I2 = 30;
    I3 = 15;
    
    Q_matrix = eye(6);
    rho = 5;
    R_matrix = rho*eye(3);
    
    inc = deg2rad(86);                  

    R_0 = 6871000;                      
    omega_o = sqrt(mu_e/R_0^3);
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                                              % Magnetic Constant
    
    
    C_io = [1, 0, 0;
            0, cos(inc), -sin(inc);
            0, sin(inc), cos(inc)];                             % Rotation from orbital to ECI
        
    C_oi = [cos(omega_o*t), -sin(omega_o*t)*cos(inc), sin(omega_o*t);
            sin(omega_o*t), cos(omega_o*t)*cos(inc), -cos(omega_o*t)*sin(inc);
            0, sin(inc), cos(inc)];                             % Rotation from ECI to perifocal
        
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];                                            % Rotation to LVLH
    C_bo = eye(3);                                              % Not body-frame in calculation of optimal gains
        
        
    R = C_io * [R_0*cos(omega_o*t); R_0*sin(omega_o*t); 0];     % Position in inertial frame
    X = R(1);
    Y = R(2);
    Z = R(3);
    
    B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];    % Magnetic field in inertial frame
    B_o = C_bo*C_oo*C_oi*B_g;                                   % Magnetic field in body frame
    B_ox = [0, -B_o(3), B_o(2); 
            B_o(3), 0, -B_o(1);
            -B_o(2), B_o(1), 0];                                % Skew-symmetric cross matrix


    Bc_lower = -I_inv*B_ox;

    Bc =[0, 0, 0;
        0, 0, 0;
        0, 0, 0;
        Bc_lower(1,1), Bc_lower(1,2), Bc_lower(1,3);
        Bc_lower(2,1), Bc_lower(2,2), Bc_lower(2,3);
        Bc_lower(3,1), Bc_lower(3,2), Bc_lower(3,3);];
    
    thetadot = Ac*theta;                                        % Open loop
end

function [Bc, B_o] = solve_Bc(t, theta)
    G = 6.674e-11;                      
    M_e = 5.97e24;                      
    mu_e = G * M_e;                     

    I1 = 22;                            
    I2 = 30;
    I3 = 15;
    
    Q_matrix = eye(6);
    rho = 5;
    R_matrix = rho*eye(3);
    
    inc = deg2rad(86);                             

    R_0 = 6871000;                      
    omega_o = sqrt(mu_e/R_0^3);
    k_1 = (I2 - I3)/I1; 
    k_3 = (I2 - I1)/I3;
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                      
    
    
    C_io = [1, 0, 0;
            0, cos(inc), -sin(inc);
            0, sin(inc), cos(inc)];
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];
        
    Bc = zeros(length(t), 6, 3);
    B_o = zeros(length(t), 3);
    
    for i = 1:length(t)
        C_oi = [cos(omega_o*t(i)), -sin(omega_o*t(i))*cos(inc), sin(omega_o*t(i));
                sin(omega_o*t(i)), cos(omega_o*t(i))*cos(inc), -cos(omega_o*t(i))*sin(inc);
                0, sin(inc), cos(inc)];
        
        C_bo = eye(3);
        
        R = C_io * [R_0*cos(omega_o*t(i)); R_0*sin(omega_o*t(i)); 0];
        
        X = R(1);
        Y = R(2);
        Z = R(3);

        B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];
        B_o(i, :) = C_bo*C_oo*C_oi*B_g;
        B_ox = [0, -B_o(i, 3), B_o(i, 2);
                B_o(i, 3), 0, -B_o(i, 1);
                -B_o(i, 2), B_o(i, 1), 0];
        Bc_lower = -I_inv*B_ox;

        Bc(i, :, :) =[0, 0, 0;
                        0, 0, 0;
                        0, 0, 0;
                        Bc_lower(1,1), Bc_lower(1,2), Bc_lower(1,3);
                        Bc_lower(2,1), Bc_lower(2,2), Bc_lower(2,3);
                        Bc_lower(3,1), Bc_lower(3,2), Bc_lower(3,3)];
    end
end

function [thetadot] = lqr_control(t, theta, Ac, Bc, K, t_list, iter, T, impulses)
    % Find the closest time on the time grid to the inputted time from
    % ode45
    t = (iter*T/(impulses-1)) + t;
    t_diff = abs(t - t_list);
    [val, idx] = min(t_diff);
    Bc_current = squeeze(Bc(idx, :, :));
    K_current = squeeze(K(idx, :, :));
    thetadot = (Ac - (Bc_current*K_current))*theta;             % Closed-loop
end
